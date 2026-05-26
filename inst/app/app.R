library(Matrix)
library(shiny)
library(bslib)
library(shinyjs)
library(memuse)
library(dplyr)
library(foreach)


options(shiny.maxRequestSize = 1000 * 1000^2)

ui <- page_fillable(
  useShinyjs(),
  theme = bs_theme(
    bg = "#0d1117",
    fg = "#e6edf3",
    primary = "#58a6ff",
    secondary = "#21262d",
    success = "#3fb950",
    font_scale = 1.1,
    bootswatch = NULL
  ),
  
  tags$head(
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Roboto+Mono:ital,wght@0,400;0,700;1,400&display=swap');
      body { font-family: 'Roboto Mono', monospace; background-color: #0d1117; color: #e6edf3; }
      .app-title { font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 3.2rem; letter-spacing: -0.03em; background: linear-gradient(90deg, #58a6ff, #3fb950); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; margin-bottom: 0; }
      .app-subtitle { font-family: 'Roboto Mono', monospace; font-size: 1.1rem; color: #8b949e; letter-spacing: 0.05em; margin-top: 4px; overflow: visible; line-height: 1.4; }
      .step-card { background: #161b22; border: 1px solid #30363d; border-radius: 12px; padding: 28px; height: 100%; position: relative; transition: border-color 0.2s ease; }
      .step-card:hover { border-color: #58a6ff44; }
      .step-badge { font-family: 'Roboto Mono', monospace; font-size: 0.65rem; font-weight: 700; letter-spacing: 0.15em; text-transform: uppercase; color: #58a6ff; background: #58a6ff15; border: 1px solid #58a6ff33; border-radius: 20px; padding: 3px 10px; display: inline-block; margin-bottom: 12px; }
      .step-badge.step2 { color: #3fb950; background: #3fb95015; border-color: #3fb95033; }
      .step-title { font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 1.1rem; margin-bottom: 6px; color: #e6edf3; }
      .step-desc { font-family: 'Roboto Mono', monospace; font-size: 0.82rem; color: #8b949e; margin-bottom: 20px; line-height: 1.6; }
      .upload-area .form-control, .upload-area input[type=file] { background: #0d1117; border: 1px dashed #30363d; border-radius: 8px; color: #8b949e; font-family: 'Roboto Mono', monospace; font-size: 0.78rem; padding: 10px; transition: border-color 0.2s; }
      .upload-area input[type=file]:hover { border-color: #58a6ff; }
      .upload-label { font-family: 'Roboto Mono', monospace; font-size: 0.72rem; color: #8b949e; text-transform: uppercase; letter-spacing: 0.08em; margin-bottom: 5px; }
      .btn-train { background: linear-gradient(135deg, #1f6feb, #388bfd); border: none; border-radius: 8px; color: white; font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 0.88rem; letter-spacing: 0.05em; padding: 10px 24px; width: 100%; margin-top: 16px; transition: all 0.2s ease; }
      .btn-train:hover { background: linear-gradient(135deg, #388bfd, #58a6ff); box-shadow: 0 0 20px #58a6ff33; transform: translateY(-1px); }
      .btn-runtime { background: linear-gradient(135deg, #6e40c9, #8957e5); border: none; border-radius: 8px; color: white; font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 0.88rem; letter-spacing: 0.05em; padding: 10px 24px; width: 100%; margin-top: 16px; transition: all 0.2s ease; }
      .btn-runtime:hover { background: linear-gradient(135deg, #8957e5, #a371f7); box-shadow: 0 0 20px #8957e533; transform: translateY(-1px); }
      .btn-run { background: linear-gradient(135deg, #238636, #2ea043); border: none; border-radius: 8px; color: white; font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 0.88rem; letter-spacing: 0.05em; padding: 10px 24px; width: 100%; margin-top: 16px; transition: all 0.2s ease; }
      .btn-run:hover { background: linear-gradient(135deg, #2ea043, #3fb950); box-shadow: 0 0 20px #3fb95033; transform: translateY(-1px); }
      .btn-home { background: linear-gradient(135deg, #1f6feb, #388bfd); border: none; border-radius: 8px; color: white; font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 0.95rem; letter-spacing: 0.05em; padding: 12px 32px; margin: 8px; transition: all 0.2s ease; }
      .btn-home:hover { transform: translateY(-1px); box-shadow: 0 0 20px #58a6ff33; }
      .btn-home.green { background: linear-gradient(135deg, #238636, #2ea043); }
      .btn-home.green:hover { box-shadow: 0 0 20px #3fb95033; }
      .btn-train:disabled, .btn-run:disabled, .btn-runtime:disabled { opacity: 0.4; transform: none; box-shadow: none; cursor: not-allowed; }
      .param-section { background: #0d1117; border: 1px solid #21262d; border-radius: 8px; padding: 14px 16px; margin-top: 16px; }
      .param-title { font-family: 'Roboto Mono', monospace; font-size: 0.68rem; color: #58a6ff; text-transform: uppercase; letter-spacing: 0.1em; margin-bottom: 12px; }
      .form-select, .form-control { background: #161b22 !important; border: 1px solid #30363d !important; color: #e6edf3 !important; border-radius: 6px !important; font-family: 'Roboto Mono', monospace !important; font-size: 0.8rem !important; }
      .form-label { font-family: 'Roboto Mono', monospace; font-size: 0.7rem; color: #8b949e; text-transform: uppercase; letter-spacing: 0.07em; }
      .status-box { background: #0d1117; border: 1px solid #30363d; border-radius: 8px; padding: 12px 16px; font-family: 'Roboto Mono', monospace; font-size: 0.75rem; color: #8b949e; margin-top: 16px; min-height: 60px; line-height: 1.7; }
      .status-success { color: #3fb950; } .status-error { color: #f85149; } .status-running { color: #d29922; }
      .header-bar { display: flex; align-items: center; justify-content: space-between; padding: 20px 28px 28px; border-bottom: 1px solid #21262d; margin-bottom: 24px; overflow: visible; }
      .author-tag { font-family: 'Roboto Mono', monospace; font-size: 0.98rem; color: #484f58; letter-spacing: 0.05em; }
      .home-card { background: #161b22; border: 1px solid #30363d; border-radius: 12px; padding: 32px; margin-bottom: 24px; }
      .home-question { font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 1.2rem; color: #e6edf3; text-align: center; margin: 24px 0 16px; }
      .home-desc { font-family: 'Roboto Mono', monospace; font-size: 0.88rem; color: #8b949e; line-height: 1.8; margin-bottom: 20px; }
      .nav-tabs { border-bottom: 1px solid #21262d !important; }
      .nav-tabs .nav-link { color: #8b949e !important; font-family: 'Roboto Mono', monospace; font-size: 0.8rem; letter-spacing: 0.05em; border: none !important; padding: 10px 20px; }
      .nav-tabs .nav-link.active { color: #58a6ff !important; background: transparent !important; border-bottom: 2px solid #58a6ff !important; }
      .nav-tabs .nav-link:hover { color: #e6edf3 !important; }
      .shiny-input-container { margin-bottom: 12px; }
      .trait-hitsize-row { display: flex; align-items: center; justify-content: space-between; padding: 6px 0; border-bottom: 1px solid #21262d; }
      .trait-name { font-family: 'Roboto Mono', monospace; font-size: 0.75rem; color: #8b949e; flex: 1; }
      .info-icon { display: inline-block; width: 14px; height: 14px; background: #30363d; border-radius: 50%; color: #8b949e; font-size: 0.6rem; font-weight: 700; text-align: center; line-height: 14px; cursor: pointer; margin-left: 5px; position: relative; vertical-align: middle; }
      .info-icon:hover .tooltip-text { visibility: visible; opacity: 1; }
      .tooltip-text { visibility: hidden; opacity: 0; background: #21262d; border: 1px solid #30363d; color: #8b949e; font-family: 'Roboto Mono', monospace; font-size: 0.68rem; border-radius: 6px; padding: 6px 10px; position: absolute; z-index: 999; bottom: 120%; left: 50%; transform: translateX(-50%); width: 200px; line-height: 1.5; transition: opacity 0.2s; pointer-events: none; }
      .cores-hint { font-family: 'Roboto Mono', monospace; font-size: 0.72rem; color: #d29922; margin-top: 5px; line-height: 1.5; }
    "))
  ),
  
  div(class = "header-bar",
      div(
        div(class = "app-title", "SVENETICS"),
        div(class = "app-subtitle", "A Multi Trait GWAS Pipeline based on Bayesian Variable Selection")
      ),
      div(class = "author-tag", "Debarshi Chakraborty, Somak Dutta, Vivekananda Roy")
  ),
  
  div(style = "padding: 0 24px 24px;",
      navset_tab(
        id = "main_tabs",
        
        nav_panel("Home",
                  div(class = "home-card", style = "margin-top: 20px;",
                      div(class = "home-desc",
                          tags$h5(style = "color: #e6edf3; font-family: 'Roboto Mono', monospace; font-weight: 700; margin-bottom: 16px;",
                                  "Welcome to SVENETICS !!!"),
                          p("SVENETICS is a Bayesian multi-trait GWAS pipeline built on top of ",
                            tags$b("SVEN"), " (Selection of Variables with Embedded screening using Bayesian methods). ",
                            "It is designed to identify causal SNPs across multiple quantitative traits in a statistically rigorous and computationally efficient manner. This web app is supposed to be a user-friendly interface to run the SVENETICS pipeline without needing to write any code. Let;s dive into the step by step details on how to use it.",),
                          tags$h5(style = "color: #e6edf3; font-family: 'Roboto Mono', monospace; font-weight: 700; margin: 20px 0 12px;",
                                  "User guide - step by step "),
                          tags$ol(style = "color: #8b949e; line-height: 2; font-size: 0.88rem;",
                                  tags$li(tags$b(style = "color:#58a6ff;", "Input SNP Matrix — "),
                                          "Supposed to be provided by the user, currently RDS format is supported only. So please convert into a .rds file before uploading. The row names of this matrix should be the genotype names. The column names should be SNP ids"),
                                  tags$li(tags$b(style = "color:#58a6ff;", " Train GWAS Model and save it — "),
                                          "Navigate to the TRAINING tab to calibrate optimal parameters for sven. An trained object will be downloaded in your computer which you will need in the multi trait GWAS step. This is independent of the phenotype and is a one time job. "),
                                  tags$li(tags$b(style = "color:#58a6ff;", "Input Trait File - "),
                                          "Supposed to be provided by the user, should be a CSV file. First column of the file should contain genotype names. The other columns may have the measurments of the different phenotypic traits."),
                                  tags$li(tags$b(style = "color:#58a6ff;", "Run Multi Trait GWAS — "),
                                          "Navigate to the RUN GWAS tab, upload the trained object and trait file to run the full multi-trait GWAS. Results will be saved in your computer as CSV files, one for each trait.")
                          ),
                          tags$h5(style = "color: #e6edf3; font-family: 'Roboto Mono', monospace; font-weight: 700; margin: 20px 0 12px;",
                                  "Important Notes"),
                          tags$ul(style = "color: #8b949e; line-height: 2; font-size: 0.88rem;",
                                  tags$li(tags$b(style = "color:#58a6ff;", "Genotype names — "),
                                          "Ideally, the genotype names of the SNPs (rownames of the .rds file) and the genotype names (first column) of the trait file should match. The rows /genotypes where names do not match will be discarded before performing GWAS."),
                                  tags$li(tags$b(style = "color:#58a6ff;", "Missing Data — "),
                                          "The SNP Matrix (in .rds format) should not contain any missing values. All missing value should be imputed."),
                                  tags$li(tags$b(style = "color:#58a6ff;", "LD Filtering — "),
                                          "SVENETICS does not perform any LD pruning/filtering. That should be done beforehand."),
                          ),
                          p(style = "color: #8b949e; font-size: 0.88rem; line-height: 1.8;",
                            "Start by training the model on your SNP matrix in the ",
                            tags$b(style = "color:#58a6ff;", "Train"), " tab. Once training is complete, upload the trained object along with your trait file in the ",
                            tags$b(style = "color:#3fb950;", "Run GWAS"), " tab to perform the full multi-trait GWAS.")
                      ),
                      div(class = "home-question", "Do you have the parameters trained?"),
                      div(style = "text-align: center;",
                          actionButton("go_train", "No — take me to training", class = "btn-home", icon = icon("dna")),
                          actionButton("go_run", "Yes — run GWAS", class = "btn-home green", icon = icon("play"))
                      )
                  )
        ),
        
        nav_panel("Train",
                  div(style = "margin-top: 20px;",
                      div(class = "step-card",
                          div(class = "step-badge", "PARAMETER TUNING"),
                          div(class = "step-title", "Train GWAS Model"),
                          div(class = "step-desc",
                              "Upload your SNP matrix (.rds) to calibrate optimal 
                               parameters. Please remember that the genotype names of this SNP.rds file should match the genotype names(first column) of the traits.csv file ideally. The names which do not match will be dropped.The trained object will be saved to your working directory."),
                          div(class = "upload-area",
                              div(class = "upload-label", "SNP Matrix (.rds)"),
                              fileInput("snp_file", label = NULL, accept = ".rds", placeholder = "X_matrix.rds", buttonLabel = "Browse")
                          ),
                          div(class = "param-section",
                              div(class = "param-title", "⚙ TRAINING PARAMETERS"),
                              div(
                                tags$label(class = "form-label", "Number of Cores",
                                           tags$span(class = "info-icon", "i",
                                                     tags$span(class = "tooltip-text", "Number of CPU cores for parallel training. 'all' uses detectCores() - 1."))),
                                textInput("n_cores", label = NULL, value = "all", placeholder = "all or e.g. 4"),
                                uiOutput("cores_hint")
                              ),
                              layout_columns(col_widths = c(6, 6),
                                             div(
                                               tags$label(class = "form-label", "Expected Hit Size",
                                                          tags$span(class = "info-icon", "i",
                                                                    tags$span(class = "tooltip-text", "Expected number of causal SNPs. 'All' tunes for small, medium and large simultaneously."))),
                                               selectInput("hitsize", label = NULL,
                                                           choices = c("All (recommended)" = "all", "Small (~10)" = "small", "Medium (~20)" = "medium", "Large (~50)" = "large"),
                                                           selected = "all")
                                             ),
                                             div(
                                               tags$label(class = "form-label", "R²",
                                                          tags$span(class = "info-icon", "i",
                                                                    tags$span(class = "tooltip-text", "Proportion of phenotypic variance explained by genetics (heritability)."))),
                                               numericInput("R2", label = NULL, value = 0.5, min = 0.01, max = 0.99, step = 0.05)
                                             )
                              ),
                              div(
                                tags$label(class = "form-label", "MAF Threshold",
                                           tags$span(class = "info-icon", "i",
                                                     tags$span(class = "tooltip-text", "Minor Allele Frequency cutoff — SNPs below this are excluded."))),
                                numericInput("MAF_threshold", label = NULL, value = 0.05, min = 0.01, max = 0.5, step = 0.01)
                              )
                          ),
                          div(
                            tags$label(class = "form-label", "Save Location",
                                       tags$span(class = "info-icon", "i",
                                                 tags$span(class = "tooltip-text", "Choose where to save the trained object."))),
                            selectInput("save_location", label = NULL,
                                        choices = c("Downloads" = "downloads", "Desktop" = "desktop", "Working Directory" = "working"),
                                        selected = "downloads"),
                            uiOutput("resolved_dir_display")
                          ),
                          actionButton("runtime_btn", "CALCULATE RUNTIME", class = "btn-runtime", icon = icon("clock")),
                          shinyjs::hidden(actionButton("train_btn", "TRAIN.GWAS", class = "btn-train", icon = icon("dna"))),
                          div(class = "status-box", uiOutput("train_status"))
                      )
                  )
        ),
        
        nav_panel("Run GWAS",
                  div(style = "margin-top: 20px;",
                      div(class = "step-card",
                          div(class = "step-badge step2", "Multi-Trait GWAS"),
                          div(class = "step-title", "Run GWAS Pipeline with optimal hyperparameters"),
                          div(class = "step-desc",
                              "Upload your saved ", tags$code("svenetics_trained"), " object (.rds) and trait file (.csv) to run the full ",
                              "multi-trait GWAS. Results for each trait will be saved as ", tags$code("{trait}_GWAS_results.csv"), " in your working directory."),
                          div(class = "upload-area",
                              div(class = "upload-label", "Trained Object (.rds)"),
                              fileInput("trained_obj_file", label = NULL, accept = ".rds", placeholder = "svenetics_trained.rds", buttonLabel = "Browse"),
                              div(class = "upload-label", style = "margin-top:10px;", "Trait File (.csv)"),
                              fileInput("trait_file", label = NULL, accept = ".csv", placeholder = "traitfile.csv", buttonLabel = "Browse")
                          ),
                          uiOutput("trait_hitsize_selectors"),
                          uiOutput("hitsize_selector"),
                          div(class = "param-section",
                              div(class = "param-title", "ℹ Pipeline Info"),
                              uiOutput("pipeline_info")
                          ),
                          actionButton("run_btn", "RUN.GWAS", class = "btn-run", icon = icon("play")),
                          div(class = "status-box", uiOutput("run_status"))
                      )
                  )
        )
      )
  )
)

server <- function(input, output, session) {
  
  trained_obj    <- reactiveVal(NULL)
  snp_matrix     <- reactiveVal(NULL)
  trait_names    <- reactiveVal(NULL)
  max_safe_cores <- reactiveVal(NULL)
  
  
  resolve_dir <- reactive({
    choice <- input$save_location
    home <- path.expand("~")
    candidate <- switch(choice,
                        "downloads" = file.path(home, "Downloads"),
                        "desktop"   = file.path(home, "Desktop"),
                        "working"   = getwd()
    )
    base <- if(dir.exists(candidate)) candidate else getwd()
    save_dir <- file.path(base, "SVENETICS_FOLDER")
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    save_dir
  })
  
  output$resolved_dir_display <- renderUI({
    tags$div(style = "font-family:'Roboto Mono',monospace; font-size:0.72rem; color:#484f58; margin-top:5px;",
             paste0("→ ", resolve_dir()))
  })
  
  
  # create output directory in user's working directory
  #out_dir <- file.path(getwd(), "SVENETICS_RESULTS")
  #if(!dir.exists(out_dir)){
    #dir.create(out_dir)
  #}
  #save_dir <- path.expand("~/SVENETICS_ RESULTS")
  

  
  get_cores <- function(){
    val      <- trimws(input$n_cores)
    max_safe <- max_safe_cores()
    if(is.null(max_safe)) max_safe <- parallel::detectCores() - 1
    if(tolower(val) == "all") return(max(1, min(parallel::detectCores() - 1, max_safe)))
    n <- suppressWarnings(as.integer(val))
    if(is.na(n) || n < 1) stop("Please enter a valid number of cores or 'all'.")
    if(n > max_safe) stop(paste0("Exceeds safe core limit. Max recommended: ", max_safe, " core(s) based on total RAM."))
    return(min(n, parallel::detectCores() - 1))
  }
  
  observeEvent(input$go_train, { nav_select("main_tabs", "Train") })
  observeEvent(input$go_run,   { nav_select("main_tabs", "Run GWAS") })
  
  observeEvent(input$snp_file, {
    req(input$snp_file)
    tryCatch({
      X <- readRDS(input$snp_file$datapath)
      X <- as(X, "dgCMatrix")
      X@x <- as.double(X@x)
      attributes(X) <- attributes(X)[c("i", "p", "Dim", "Dimnames", "x", "factors", "class")]
      snp_matrix(X)
      shinyjs::hide("train_btn")
      mat_size_mb  <- as.numeric(object.size(X)) / (1000^2)
      total_ram_gb <- as.numeric(memuse::Sys.meminfo()$totalram)/(1000^3)
      total_ram_mb <- total_ram_gb * 1000
      #total_ram_mb <- as.numeric(memuse::Sys.meminfo()$totalram, units = "MiB")
      safe_cores   <- max(1, min(floor(total_ram_mb / mat_size_mb), parallel::detectCores() - 1))
      safe_cores <- min(100,safe_cores)
      max_safe_cores(safe_cores)
      output$cores_hint <- renderUI({
        tags$div(class = "cores-hint",
                 paste0("⚠ Max recommended cores: ", safe_cores), tags$br(),
                 paste0("Matrix: ", round(mat_size_mb, 1), " MB  |  Total RAM: ", round(total_ram_gb, 0), " GB"))
        #paste0("Matrix: ", round(mat_size_mb, 1), " MB  |  Total RAM: ", round(total_ram_mb, 0), " MB"))
      })
    }, error = function(e) {
      output$train_status <- renderUI({
        tags$span(class = "status-error", paste("✖ Could not load SNP matrix:", e$message))
      })
    })
  })
  
  observeEvent(input$trait_file, {
    if(is.null(input$trait_file)) return()
    tryCatch({
      tf <- read.csv(input$trait_file$datapath, row.names = NULL)
      trait_names(colnames(tf)[-1])
    }, error = function(e) {
      output$run_status <- renderUI({
        tags$span(class = "status-error", paste("✖ Could not read trait file:", e$message))
      })
    })
  })
  
  output$trait_hitsize_selectors <- renderUI({
    tn  <- trait_names()
    obj <- trained_obj()
    if(is.null(tn) || is.null(obj)) return(NULL)
    is_all <- is.list(obj$optimal_params) && all(c("small","medium","large") %in% names(obj$optimal_params))
    if(!is_all) return(NULL)
    div(class = "param-section", style = "margin-top: 16px;",
        div(class = "param-title", "⚙ EXPECTED HIT SIZE PER TRAIT"),
        tagList(lapply(seq_along(tn), function(i){
          div(class = "trait-hitsize-row",
              div(class = "trait-name", tn[i]),
              div(style = "width: 160px;",
                  selectInput(inputId = paste0("hitsize_trait_", i), label = NULL,
                              choices = c("Small (~10)" = "small", "Medium (~20)" = "medium", "Large (~50)" = "large"),
                              selected = "medium")))
        }))
    )
  })
  
  output$hitsize_selector <- renderUI({
    tn  <- trait_names()
    obj <- trained_obj()
    if(is.null(obj)) return(NULL)
    if(!is.null(tn)) return(NULL)
    is_all <- is.list(obj$optimal_params) && all(c("small","medium","large") %in% names(obj$optimal_params))
    if(!is_all) return(NULL)
    div(class = "param-section", style = "margin-bottom: 16px;",
        div(class = "param-title", "⚙ SELECT HIT SIZE FOR GWAS"),
        selectInput("run_hitsize", label = NULL,
                    choices = c("Small (~10)" = "small", "Medium (~20)" = "medium", "Large (~50)" = "large"),
                    selected = "medium"))
  })
  
  observeEvent(input$runtime_btn, {
    req(snp_matrix())
    shinyjs::disable("runtime_btn")
    shinyjs::html("runtime_btn", "<i class='fa fa-clock'></i> Calculating training time...")
    on.exit({ shinyjs::enable("runtime_btn"); shinyjs::html("runtime_btn", "<i class='fa fa-clock'></i> CALCULATE RUNTIME") })
    tryCatch({
      X       <- snp_matrix()
      n_cores <- get_cores()
      t_single <- as.numeric(calc.runtime(X), units = "mins")
      if(input$hitsize == "all"){
        t_total   <- round(t_single * 3600 * 3 / n_cores, 2)
        size_note <- " (tuning all 3 hit sizes)"
      } else {
        t_total   <- round(t_single * 3600 / n_cores, 2)
        size_note <- ""
      }
      output$train_status <- renderUI({
        tags$span(class = "status-success",
                  paste0("⏱ Single SVEN runtime: ", round(t_single, 3), " mins"), tags$br(),
                  paste0("📊 Estimated training time: ", t_total, " mins on ", n_cores, " core(s)", size_note), tags$br(),
                  "Click the TRAIN.GWAS button located above to start.")
      })
      shinyjs::show("train_btn")
    }, error = function(e) {
      output$train_status <- renderUI({ tags$span(class = "status-error", paste("✖ Error:", e$message)) })
    })
  })
  
  observeEvent(input$train_btn, {
    req(snp_matrix())
    shinyjs::disable("train_btn")
    shinyjs::html("train_btn", "<i class='fa fa-dna'></i> Training started...")
    on.exit({ shinyjs::enable("train_btn"); shinyjs::html("train_btn", "<i class='fa fa-dna'></i> TRAIN.GWAS") })
    tryCatch({
      X       <- snp_matrix()
      n_cores <- get_cores()
      obj <- isolate({
        parameter_selection(X = X, R2 = input$R2, betamax = 1, n.cores = n_cores,
                            hitsize = input$hitsize, MAF_threshold = input$MAF_threshold)
      })
      trained_obj(obj)
      dir.create(resolve_dir(), showWarnings = FALSE, recursive = TRUE)
      save_path <- file.path(resolve_dir(), "svenetics_trained.rds")
      saveRDS(obj, file = save_path)
      if(input$hitsize == "all"){
        param_text <- tagList(
          "Small  → λ = ", round(obj$optimal_params$small["lambda"], 4), " | w = ", round(obj$optimal_params$small["w"], 4), tags$br(),
          "Medium → λ = ", round(obj$optimal_params$medium["lambda"], 4), " | w = ", round(obj$optimal_params$medium["w"], 4), tags$br(),
          "Large  → λ = ", round(obj$optimal_params$large["lambda"], 4), " | w = ", round(obj$optimal_params$large["w"], 4))
      } else {
        param_text <- tagList("Optimal λ = ", round(obj$optimal_params["lambda"], 4), " | w = ", round(obj$optimal_params["w"], 4))
      }
      output$train_status <- renderUI({
        tags$span(class = "status-success", "✔ Training complete!", tags$br(), param_text, tags$br(),
                  paste0("Saved → ", save_path))
      })
    }, error = function(e) {
      output$train_status <- renderUI({ tags$span(class = "status-error", paste("✖ Error:", e$message)) })
    })
  })
  
  observeEvent(input$trained_obj_file, {
    if(is.null(input$trained_obj_file)) return()
    tryCatch({
      obj <- readRDS(input$trained_obj_file$datapath)
      trained_obj(obj)
    }, error = function(e) {
      output$run_status <- renderUI({ tags$span(class = "status-error", paste("✖ Could not load object:", e$message)) })
    })
  })
  
  output$pipeline_info <- renderUI({
    obj <- trained_obj()
    if(is.null(obj)){
      return(tags$span(style = "color:#484f58; font-family:'Roboto Mono',monospace; font-size:0.75rem;", "Upload a trained object to see details."))
    }
    is_all <- is.list(obj$optimal_params) && all(c("small","medium","large") %in% names(obj$optimal_params))
    tags$div(style = "font-family:'Roboto Mono',monospace; font-size:0.75rem; color:#8b949e; line-height:2;",
             if(is_all){ tagList(tags$span(style = "color:#58a6ff;", "Tuned for "), "small, medium, large", tags$br())
             } else { tagList(tags$span(style = "color:#58a6ff;", "λ "), round(obj$optimal_params["lambda"], 4), tags$br(),
                              tags$span(style = "color:#58a6ff;", "w "), round(obj$optimal_params["w"], 4), tags$br()) },
             tags$span(style = "color:#58a6ff;", "SNPs (cleaned) "), ncol(obj$X), tags$br(),
             tags$span(style = "color:#58a6ff;", "Cores "), obj$n.cores, tags$br(),
             tags$span(style = "color:#58a6ff;", "Traits "), if(!is.null(trait_names())) length(trait_names()) else "supplied at run time"
    )
  })
  
  observeEvent(input$run_btn, {
    req(trained_obj(), input$trait_file)
    #------------------------------------------
    shinyjs::disable("run_btn") 
    shinyjs::html("run_btn", "<i class='fa fa-dna'></i> ⏳ Running GWAS pipeline across all traits — please wait...") 
    on.exit({ shinyjs::enable("run_btn"); shinyjs::html("run_btn", "<i class='fa fa-dna'></i> RUN.GWAS") }) 
    #------------------------------------------
    output$run_status <- renderUI({
      tags$span(class = "status-running", "⏳ Running GWAS pipeline across all traits — please wait...")
    })
    tryCatch({
      obj       <- trained_obj()
      tn        <- trait_names()
      traitfile <- read.csv(input$trait_file$datapath, row.names = NULL)
      is_all <- is.list(obj$optimal_params) && all(c("small","medium","large") %in% names(obj$optimal_params))
      if(is_all && !is.null(tn)){
        hitsizes <- sapply(seq_along(tn), function(i){
          val <- input[[paste0("hitsize_trait_", i)]]
          if(is.null(val)) "medium" else val
        })
      } else {
        hitsizes <- NULL
      }
      dir.create(resolve_dir(), showWarnings = FALSE, recursive = TRUE)
      result   <- svenetics_pipeline(obj, traitfile, hitsizes = hitsizes, save_dir = resolve_dir())
      n_traits <- length(result)
      output$run_status <- renderUI({
        tags$span(class = "status-success",
                  paste0("✔ GWAS complete! ", n_traits, " trait(s) processed."), tags$br(),
                  paste0("Results saved to: ",resolve_dir()))
      })
    }, error = function(e) {
      output$run_status <- renderUI({
        tags$span(class = "status-error", paste("✖ Error:", e$message))
      })
    })
  })
}

shinyApp(ui = ui, server = server)
