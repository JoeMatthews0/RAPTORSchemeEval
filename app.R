library(shiny)
library(DT)
library(ggplot2)
library(MASS)

# ── Helpers ────────────────────────────────────────────────────────────────────

read_upload <- function(path, filename = "") {

  # ── CSV branch ───────────────────────────────────────────────────────────────
  is_csv <- grepl("\\.csv$", filename, ignore.case = TRUE)
  if (!is_csv) {
    # Detect from content: if the first non-empty line has commas, treat as CSV
    first_line <- trimws(readLines(path, n = 1L, warn = FALSE))
    is_csv <- nchar(first_line) > 0 && grepl(",", first_line)
  }

  if (is_csv) {
    df <- tryCatch(
      read.csv(path, header = TRUE, stringsAsFactors = FALSE,
               check.names = TRUE),
      error = function(e) data.frame()
    )
    return(df)
  }

  # ── Whitespace-delimited branch (original .txt logic) ────────────────────────
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2) {
    return(data.frame())
  }

  # --- Parse header line ---
  header_raw <- trimws(lines[1])
  header_tokens <- strsplit(header_raw, "\\s+")[[1]]
  # Strip leading '#' from first token (e.g. "#Site" -> "Site")
  header_tokens[1] <- sub("^#", "", header_tokens[1])
  # Strip surrounding quotes (e.g. '"Sp40"' -> 'Sp40')
  header_tokens <- gsub('^"|"$', "", header_tokens)

  # --- Count columns in data rows ---
  data_lines <- lines[2:min(5, length(lines))]
  data_ncols <- max(lengths(strsplit(trimws(data_lines), "\\s+")), na.rm = TRUE)

  # --- Merge header tokens until count matches data ---
  # Merge tokens whose value looks like a secondary word
  # (e.g. "No.", "type", "Class") with the preceding token.
  secondary <- c(
    "No.", "no.", "type", "Type", "Class", "class",
    "Limit", "limit", "Road", "road"
  )
  toks <- header_tokens
  while (length(toks) > data_ncols) {
    idx <- which(toks %in% secondary)
    if (length(idx) > 0 && idx[1] > 1) {
      i <- idx[1]
      toks[i - 1] <- paste(toks[i - 1], toks[i], sep = "_")
      toks <- toks[-i]
    } else {
      # Fallback: merge last two tokens
      n <- length(toks)
      toks[n - 1] <- paste(toks[n - 1], toks[n], sep = "_")
      toks <- toks[-n]
    }
  }

  # --- Handle row-names case: data has 1 more token than header ---
  if (data_ncols == length(toks) + 1L) {
    toks <- c("SiteID", toks)
  }

  col_names <- make.names(toks, unique = TRUE)

  # --- Read data rows with cleaned column names ---
  df <- tryCatch(
    {
      read.table(
        text = paste(lines[-1], collapse = "\n"),
        header = FALSE,
        col.names = col_names,
        fill = TRUE,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      # Fall back to standard read; restore row-names if needed
      d <- read.table(path,
        header = TRUE, fill = TRUE,
        comment.char = "", check.names = TRUE,
        stringsAsFactors = FALSE
      )
      rn <- rownames(d)
      if (!all(rn == as.character(seq_len(nrow(d))))) {
        d <- cbind(SiteID = rn, d, stringsAsFactors = FALSE)
        rownames(d) <- NULL
      }
      d
    }
  )
  df
}

# ── Fully Bayesian MCMC ────────────────────────────────────────────────────────
#
# Model (Fawcett & Thorpe 2013):
#   Reference sites  : y_ref ~ NB(γ, μ_ref)   μ_ref = exp(X_ref β)
#   Treated sites    : y_before ~ Poisson(mj)
#                      mj ~ Gamma(γ, γ/μj)     μj = exp(X_trt β)
#   Posterior of mj  : Gamma(γ + y_before, γ/μj + 1)   [Gibbs step]
#   Priors           : βi ~ N(0, 10²),  ρ = log(κ) ~ N(0, 10²),  γ = exp(-ρ)
#   Inference for β,ρ: component-wise random-walk Metropolis-Hastings
#                      using NB log-likelihood on reference data

run_fb_mcmc <- function(X_ref, y_ref,
                        X_trt, y_before,
                        n_iter,
                        progress_fn = NULL) {
  n_beta <- ncol(X_ref)
  n_trt <- nrow(X_trt)

  # ── Starting values from NB GLM ────────────────────────────────────────────
  beta_cur <- rep(0, n_beta)
  rho_cur <- 0 # log(kappa); kappa = 1/gamma

  tryCatch(
    {
      df_fit <- as.data.frame(X_ref[, -1, drop = FALSE]) # drop intercept
      df_fit$y <- y_ref
      fit <- glm.nb(y ~ ., data = df_fit)
      beta_cur <- as.numeric(coef(fit))
      rho_cur <- log(1 / fit$theta)
    },
    error = function(e) NULL
  )

  # ── Log-posterior (reference data likelihood + priors) ─────────────────────
  log_post <- function(beta, rho) {
    gamma <- exp(-rho)
    mu <- as.numeric(exp(X_ref %*% beta))
    if (any(!is.finite(mu)) || gamma <= 0) {
      return(-Inf)
    }
    ll <- sum(dnbinom(y_ref, size = gamma, mu = mu, log = TRUE))
    if (!is.finite(ll)) {
      return(-Inf)
    }
    lp <- sum(dnorm(beta, 0, 10, log = TRUE)) + dnorm(rho, 0, 10, log = TRUE)
    ll + lp
  }

  # ── Storage ────────────────────────────────────────────────────────────────
  m_samples <- matrix(NA_real_, nrow = n_iter, ncol = n_trt)

  # ── Adaptive MH tuning ─────────────────────────────────────────────────────
  sigma_beta <- rep(0.1, n_beta)
  sigma_rho <- 0.2
  acc_beta <- integer(n_beta)
  acc_rho <- 0L
  tune_every <- 100L
  burn_in <- max(1L, as.integer(n_iter * 0.2))

  lp_cur <- log_post(beta_cur, rho_cur)

  # ── MCMC loop ──────────────────────────────────────────────────────────────
  for (iter in seq_len(n_iter)) {
    if (!is.null(progress_fn) && iter %% 200L == 0L) {
      progress_fn(iter / n_iter)
    }

    # Component-wise MH for each beta
    for (k in seq_len(n_beta)) {
      bp <- beta_cur
      bp[k] <- beta_cur[k] + rnorm(1L, 0, sigma_beta[k])
      lp_prop <- log_post(bp, rho_cur)
      if (is.finite(lp_prop) && log(runif(1L)) < lp_prop - lp_cur) {
        beta_cur <- bp
        lp_cur <- lp_prop
        if (iter <= burn_in) acc_beta[k] <- acc_beta[k] + 1L
      }
    }

    # MH for rho
    rp <- rho_cur + rnorm(1L, 0, sigma_rho)
    lp_prop <- log_post(beta_cur, rp)
    if (is.finite(lp_prop) && log(runif(1L)) < lp_prop - lp_cur) {
      rho_cur <- rp
      lp_cur <- lp_prop
      if (iter <= burn_in) acc_rho <- acc_rho + 1L
    }

    # Adaptive tuning during burn-in
    if (iter <= burn_in && iter %% tune_every == 0L) {
      for (k in seq_len(n_beta)) {
        rate <- acc_beta[k] / tune_every
        if (rate > 0.44) {
          sigma_beta[k] <- sigma_beta[k] * 1.3
        } else if (rate < 0.20) sigma_beta[k] <- sigma_beta[k] * 0.7
        acc_beta[k] <- 0L
      }
      rate_rho <- acc_rho / tune_every
      if (rate_rho > 0.44) {
        sigma_rho <- sigma_rho * 1.3
      } else if (rate_rho < 0.20) sigma_rho <- sigma_rho * 0.7
      acc_rho <- 0L
    }

    # Gibbs sample mj for each treated site
    gamma_cur <- exp(-rho_cur)
    mu_trt <- as.numeric(exp(X_trt %*% beta_cur))
    shape_j <- gamma_cur + y_before
    rate_j <- gamma_cur / mu_trt + 1
    m_samples[iter, ] <- rgamma(n_trt, shape = shape_j, rate = rate_j)
  }

  # Discard burn-in, return posterior means
  keep <- seq(burn_in + 1L, n_iter)
  m_post <- m_samples[keep, , drop = FALSE]
  colMeans(m_post)
}

# ── UI ─────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; margin-top: 0; }
    .raptor-navbar {
      background-color: #2c3e50;
      border-bottom: 3px solid #e1202b;
      padding: 0 20px;
      margin-bottom: 20px;
      min-height: 50px;
      display: flex;
      align-items: center;
    }
    .raptor-navbar .raptor-brand {
      color: #ffffff;
      font-size: 18px;
      font-weight: bold;
      letter-spacing: 0.5px;
      margin: 0;
      padding: 13px 0;
    }
    .raptor-navbar .raptor-brand span {
      color: #e1202b;
    }
    .sidebar-panel-custom { background: #f5f5f5; padding: 15px; border-radius: 4px; }
    h4 { color: #2c3e50; border-bottom: 2px solid #e1202b; padding-bottom: 4px; margin-top: 18px; }
    .btn-run { width: 100%; margin-top: 12px; }
    .btn-primary { background-color: #e1202b; border-color: #c0392b; }
    .btn-primary:hover, .btn-primary:focus {
      background-color: #c0392b; border-color: #a93226; }
  "))),
  div(class = "raptor-navbar",
    p(class = "raptor-brand",
      tags$span("RAPTOR"), " Scheme Evaluation")
  ),
  sidebarLayout(
    sidebarPanel(
      div(
        class = "sidebar-panel-custom",
        div(style = "margin-bottom: 8px;",
          tags$b("Quick start"),
          br(),
          actionButton("use_example_btn", "Load example data",
                       icon  = icon("database"),
                       class = "btn-default btn-sm",
                       style = "width:100%; margin-top:4px;")
        ),
        tags$hr(style = "margin: 10px 0;"),
        h4("Reference Data"),
        fileInput("ref_file", "Upload reference file (.txt / .csv)",
          accept = c(".txt", ".csv"),
          width  = "100%"
        ),
        uiOutput("ref_id_ui"),
        uiOutput("ref_count_ui"),
        uiOutput("ref_cov_ui"),
        h4("Treated Data"),
        fileInput("trt_file", "Upload treated file (.txt / .csv)",
          accept = c(".txt", ".csv"),
          width  = "100%"
        ),
        uiOutput("trt_id_ui"),
        uiOutput("trt_before_ui"),
        uiOutput("trt_after_ui"),
        helpText("Tip: select multiple columns to sum them (e.g. Fatal + Serious + Slight)."),
        uiOutput("cov_match_ui"),
        h4("MCMC Settings"),
        numericInput("n_iter", "Posterior samples",
          value = 10000, min = 2000, max = 500000, step = 1000
        ),
        actionButton("run_btn", "Run Analysis",
          class = "btn-primary btn-run",
          icon  = icon("play")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Data Preview",
          br(),
          h5("Reference Data"),
          DTOutput("ref_preview"),
          br(),
          h5("Treated Data"),
          DTOutput("trt_preview")
        ),
        tabPanel(
          "Results Table",
          br(),
          uiOutput("results_summary"),
          br(),
          DTOutput("results_table")
        ),
        tabPanel(
          "Results Plot",
          br(),
          h5("Summary: RTM vs Treatment effect by site"),
          plotOutput("results_plot", height = "400px"),
          br(),
          h5("Individual site breakdown"),
          uiOutput("site_plots_ui")
        )
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Example-data toggle ───────────────────────────────────────────────────
  using_example <- reactiveVal(FALSE)

  observeEvent(input$use_example_btn, {
    using_example(TRUE)
  })
  # Reset flag if the user uploads their own files
  observeEvent(input$ref_file, { using_example(FALSE) })
  observeEvent(input$trt_file, { using_example(FALSE) })

  # ── Reactive data ─────────────────────────────────────────────────────────

  ref_data <- reactive({
    if (using_example()) {
      read_upload("ExampleReference.csv", "ExampleReference.csv")
    } else {
      req(input$ref_file)
      read_upload(input$ref_file$datapath, input$ref_file$name)
    }
  })

  trt_data <- reactive({
    if (using_example()) {
      read_upload("ExampleTreated.csv", "ExampleTreated.csv")
    } else {
      req(input$trt_file)
      read_upload(input$trt_file$datapath, input$trt_file$name)
    }
  })

  # ── Data previews ─────────────────────────────────────────────────────────

  output$ref_preview <- renderDT({
    req(ref_data())
    datatable(head(ref_data(), 10),
      options = list(scrollX = TRUE, pageLength = 5, dom = "t"),
      rownames = FALSE
    )
  })

  output$trt_preview <- renderDT({
    req(trt_data())
    datatable(head(trt_data(), 10),
      options = list(scrollX = TRUE, pageLength = 5, dom = "t"),
      rownames = FALSE
    )
  })

  # ── Reference column selectors ────────────────────────────────────────────

  output$ref_id_ui <- renderUI({
    req(ref_data())
    cols <- colnames(ref_data())
    selectInput("ref_id_col", "Site ID column",
      choices = cols, selected = cols[1]
    )
  })

  output$ref_count_ui <- renderUI({
    req(ref_data())
    cols <- colnames(ref_data())
    selectInput("ref_count_col", "Collision count column",
      choices = cols, selected = cols[length(cols)]
    )
  })

  ref_cov_choices <- reactive({
    req(ref_data(), input$ref_id_col, input$ref_count_col)
    setdiff(colnames(ref_data()), c(input$ref_id_col, input$ref_count_col))
  })

  output$ref_cov_ui <- renderUI({
    req(ref_cov_choices())
    ch <- ref_cov_choices()
    checkboxGroupInput("ref_cov_cols", "Covariates for SPF",
      choices = ch, selected = ch
    )
  })

  # ── Treated column selectors ──────────────────────────────────────────────

  output$trt_id_ui <- renderUI({
    req(trt_data())
    cols <- colnames(trt_data())
    selectInput("trt_id_col", "Site ID column",
      choices = cols, selected = cols[1]
    )
  })

  output$trt_before_ui <- renderUI({
    req(trt_data(), input$trt_id_col)
    cols <- setdiff(colnames(trt_data()), input$trt_id_col)
    default <- if (using_example() && "Before" %in% cols) {
      "Before"
    } else {
      cols[min(2, length(cols))]
    }
    selectInput("trt_before_col", "Before count column(s)",
      choices  = cols,
      selected = default,
      multiple = TRUE
    )
  })

  output$trt_after_ui <- renderUI({
    req(trt_data(), input$trt_id_col)
    cols <- setdiff(colnames(trt_data()), input$trt_id_col)
    default <- if (using_example() && "After" %in% cols) {
      "After"
    } else {
      cols[min(5, length(cols))]
    }
    selectInput("trt_after_col", "After count column(s)",
      choices  = cols,
      selected = default,
      multiple = TRUE
    )
  })

  # ── Covariate matching (ref → treated) ────────────────────────────────────
  # For each selected reference covariate, show a dropdown to pick the
  # corresponding column in the treated dataset.

  output$cov_match_ui <- renderUI({
    req(input$ref_cov_cols, trt_data())
    sel <- input$ref_cov_cols
    if (length(sel) == 0) {
      return(NULL)
    }
    trt_cols <- setdiff(
      colnames(trt_data()),
      c(
        input$trt_id_col,
        input$trt_before_col,
        input$trt_after_col
      )
    )
    tagList(
      tags$b("Match covariates to treated data columns:"),
      lapply(sel, function(cov) {
        default <- if (cov %in% trt_cols) cov else trt_cols[1]
        selectInput(
          inputId  = paste0("trt_cov__", make.names(cov)),
          label    = cov,
          choices  = trt_cols,
          selected = default
        )
      })
    )
  })

  # ── Run analysis ──────────────────────────────────────────────────────────

  results <- eventReactive(input$run_btn, {
    req(
      ref_data(), trt_data(),
      input$ref_id_col, input$ref_count_col, input$ref_cov_cols,
      input$trt_id_col, input$trt_before_col, input$trt_after_col
    )

    ref <- ref_data()
    trt <- trt_data()

    # Collect selected reference covariates and their treated-data matches
    ref_cov <- input$ref_cov_cols
    trt_cov <- vapply(ref_cov, function(cov) {
      v <- input[[paste0("trt_cov__", make.names(cov))]]
      if (is.null(v)) cov else v
    }, character(1))

    # Build design matrices (with intercept)
    to_num <- function(x) suppressWarnings(as.numeric(x))

    X_ref_raw <- sapply(ref_cov, function(c) to_num(ref[[c]]))
    X_trt_raw <- sapply(trt_cov, function(c) to_num(trt[[c]]))

    X_ref <- cbind(1, X_ref_raw)
    X_trt <- cbind(1, X_trt_raw)

    col_sum <- function(df, cols) {
      if (length(cols) == 1) {
        return(to_num(df[[cols]]))
      }
      rowSums(sapply(cols, function(c) to_num(df[[c]])), na.rm = FALSE)
    }
    y_ref <- to_num(ref[[input$ref_count_col]])
    y_before <- col_sum(trt, input$trt_before_col)
    y_after <- col_sum(trt, input$trt_after_col)

    # Drop rows with NA
    ref_ok <- complete.cases(X_ref, y_ref)
    trt_ok <- complete.cases(X_trt, y_before, y_after)

    if (sum(ref_ok) < 5) {
      stop("Fewer than 5 complete reference-site rows — please check column selections.")
    }
    if (sum(trt_ok) < 1) {
      stop("No complete treated-site rows — please check column selections.")
    }

    X_ref <- X_ref[ref_ok, , drop = FALSE]
    y_ref <- y_ref[ref_ok]
    X_trt <- X_trt[trt_ok, , drop = FALSE]
    y_before <- y_before[trt_ok]
    y_after <- y_after[trt_ok]
    site_ids <- trt[[input$trt_id_col]][trt_ok]

    # Compute initial SPF predictions for diagnostics (uses GLM starting values)
    tryCatch(
      {
        df_fit_diag <- as.data.frame(X_ref[, -1, drop = FALSE])
        df_fit_diag$y <- y_ref
        fit_diag <- glm.nb(y ~ ., data = df_fit_diag)
        mu_trt_diag <- as.numeric(exp(X_trt %*% coef(fit_diag)))
        message("=== SPF diagnostic ===")
        message(sprintf(
          "  y_ref   range : %.1f – %.1f  (mean %.1f)",
          min(y_ref), max(y_ref), mean(y_ref)
        ))
        message(sprintf(
          "  y_before range: %.1f – %.1f  (mean %.1f)",
          min(y_before), max(y_before), mean(y_before)
        ))
        message(sprintf(
          "  mu_trt  range : %.1f – %.1f  (mean %.1f)",
          min(mu_trt_diag), max(mu_trt_diag), mean(mu_trt_diag)
        ))
        message(sprintf(
          "  sites where mu_trt > y_before: %d / %d",
          sum(mu_trt_diag > y_before), length(y_before)
        ))
      },
      error = function(e) NULL
    )

    withProgress(message = "Running Fully Bayesian MCMC…", value = 0, {
      m_mean <- run_fb_mcmc(
        X_ref = X_ref,
        y_ref = y_ref,
        X_trt = X_trt,
        y_before = y_before,
        n_iter = input$n_iter,
        progress_fn = function(p) setProgress(p)
      )
    })

    data.frame(
      Site             = site_ids,
      Before           = y_before,
      After            = y_after,
      RTM_estimate     = round(m_mean, 2),
      Raw_change       = y_after - y_before,
      RTM_change       = round(m_mean - y_before, 2),
      Treatment_change = round(y_after - m_mean, 2),
      stringsAsFactors = FALSE
    )
  })

  # ── Results summary banner ────────────────────────────────────────────────

  output$results_summary <- renderUI({
    req(results())
    df <- results()
    total_raw <- sum(df$Raw_change)
    total_rtm <- sum(df$RTM_change)
    total_trt <- sum(df$Treatment_change)
    div(
      style = "display:flex; gap:20px;",
      div(
        style = "background:#eef2ff; padding:12px 20px; border-radius:4px; flex:1; text-align:center;",
        div(style = "font-size:1.4em; font-weight:bold;", total_raw),
        div("Total raw change")
      ),
      div(
        style = "background:#fff7e6; padding:12px 20px; border-radius:4px; flex:1; text-align:center;",
        div(style = "font-size:1.4em; font-weight:bold;", round(total_rtm, 1)),
        div("Change due to RTM")
      ),
      div(
        style = "background:#edfff0; padding:12px 20px; border-radius:4px; flex:1; text-align:center;",
        div(style = "font-size:1.4em; font-weight:bold;", round(total_trt, 1)),
        div("Change due to treatment")
      )
    )
  })

  # ── Results table ─────────────────────────────────────────────────────────

  output$results_table <- renderDT({
    req(results())
    df <- results()
    names(df) <- c(
      "Site", "Before", "After", "RTM Estimate (post. mean)",
      "Raw Change", "RTM Change", "Treatment Change"
    )
    datatable(df,
      rownames = FALSE,
      options  = list(pageLength = 20, scrollX = TRUE)
    ) |>
      formatStyle("RTM Change",
        color = styleInterval(0, c("green", "red"))
      ) |>
      formatStyle("Treatment Change",
        color = styleInterval(0, c("green", "red"))
      )
  })

  # ── Results plot ──────────────────────────────────────────────────────────

  output$results_plot <- renderPlot({
    req(results())
    df <- results()

    plot_df <- rbind(
      data.frame(Site = as.character(df$Site),
                 Component = "RTM Change",
                 Value = df$RTM_change),
      data.frame(Site = as.character(df$Site),
                 Component = "Treatment Change",
                 Value = df$Treatment_change)
    )
    plot_df$Site      <- factor(plot_df$Site, levels = as.character(df$Site))
    plot_df$Component <- factor(plot_df$Component,
                                levels = c("RTM Change", "Treatment Change"))

    ggplot(plot_df, aes(x = Site, y = Value, fill = Component)) +
      geom_col(position = "stack", width = 0.7) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
      scale_fill_manual(values = c("RTM Change"       = "#5B9BD5",
                                   "Treatment Change" = "#ED7D31")) +
      labs(title = "Summary: RTM vs Treatment effect by site",
           x = "Site", y = "Change in collisions", fill = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x   = element_text(angle = 90, hjust = 1,
                                         vjust = 0.5, size = 7),
            legend.position = "bottom",
            plot.title      = element_text(face = "bold", size = 13))
  })

  # ── Individual site breakdown (faceted) ───────────────────────────────────

  # Dynamic height: 4 columns, ~200 px per row, min 250 px
  output$site_plots_ui <- renderUI({
    req(results())
    n_sites  <- nrow(results())
    n_cols   <- min(4L, n_sites)
    n_rows   <- ceiling(n_sites / n_cols)
    px_height <- max(250L, n_rows * 200L)
    plotOutput("site_facet_plot", height = paste0(px_height, "px"))
  })

  output$site_facet_plot <- renderPlot({
    req(results())
    df <- results()

    bar_colours <- c("Overall" = "#2c3e50",
                     "RTM"     = "#5B9BD5",
                     "Treatment" = "#ED7D31")

    plot_df <- rbind(
      data.frame(Site   = as.character(df$Site),
                 Metric = "Overall",
                 Value  = df$Raw_change),
      data.frame(Site   = as.character(df$Site),
                 Metric = "RTM",
                 Value  = df$RTM_change),
      data.frame(Site   = as.character(df$Site),
                 Metric = "Treatment",
                 Value  = df$Treatment_change)
    )
    plot_df$Site   <- factor(plot_df$Site,   levels = as.character(df$Site))
    plot_df$Metric <- factor(plot_df$Metric,
                             levels = c("Overall", "RTM", "Treatment"))

    # Label y position: midpoint inside bar; nudge tiny bars just outside
    plot_df$label_y <- ifelse(
      abs(plot_df$Value) < 0.5,
      ifelse(plot_df$Value >= 0, 0.3, -0.3),
      plot_df$Value / 2
    )

    n_cols <- min(4L, nrow(df))

    ggplot(plot_df, aes(x = Metric, y = Value, fill = Metric)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
      geom_text(aes(y = label_y, label = Metric),
                colour = "white", fontface = "bold", size = 3) +
      facet_wrap(~ Site, ncol = n_cols, scales = "fixed") +
      scale_fill_manual(values = bar_colours) +
      labs(x = NULL, y = "Change in collisions") +
      theme_bw(base_size = 10) +
      theme(axis.text.x       = element_blank(),
            axis.ticks.x      = element_blank(),
            strip.background  = element_rect(fill = "#2c3e50"),
            strip.text        = element_text(colour = "white", face = "bold",
                                             size = 9),
            panel.spacing     = unit(0.8, "lines"))
  })
}

shinyApp(ui = ui, server = server)
