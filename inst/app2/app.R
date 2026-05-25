library(shiny)
library(bslib)
library(Matrix)
library(shinyjs)


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
      .app-subtitle { font-family: 'Roboto Mono', monospace; font-size: 1.1rem; color: #8b949e; letter-spacing: 0.05em; margin-top: 4px; line-height: 1.4; }
      .header-bar { display: flex; align-items: center; justify-content: space-between; padding: 20px 28px 28px; border-bottom: 1px solid #21262d; margin-bottom: 24px; }
      .author-tag { font-family: 'Roboto Mono', monospace; font-size: 0.98rem; color: #484f58; letter-spacing: 0.05em; }
      .step-card { background: #161b22; border: 1px solid #30363d; border-radius: 12px; padding: 28px; position: relative; transition: border-color 0.2s ease; }
      .step-card:hover { border-color: #58a6ff44; }
      .step-badge { font-family: 'Roboto Mono', monospace; font-size: 0.65rem; font-weight: 700; letter-spacing: 0.15em; text-transform: uppercase; color: #58a6ff; background: #58a6ff15; border: 1px solid #58a6ff33; border-radius: 20px; padding: 3px 10px; display: inline-block; margin-bottom: 12px; }
      .step-title { font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 1.1rem; margin-bottom: 6px; color: #e6edf3; }
      .step-desc { font-family: 'Roboto Mono', monospace; font-size: 0.82rem; color: #8b949e; margin-bottom: 20px; line-height: 1.6; }
      .upload-area .form-control, .upload-area input[type=file] { background: #0d1117; border: 1px dashed #30363d; border-radius: 8px; color: #8b949e; font-family: 'Roboto Mono', monospace; font-size: 0.78rem; padding: 10px; transition: border-color 0.2s; }
      .upload-area input[type=file]:hover { border-color: #58a6ff; }
      .upload-label { font-family: 'Roboto Mono', monospace; font-size: 0.72rem; color: #8b949e; text-transform: uppercase; letter-spacing: 0.08em; margin-bottom: 5px; }
      .param-section { background: #0d1117; border: 1px solid #21262d; border-radius: 8px; padding: 14px 16px; margin-top: 16px; }
      .param-title { font-family: 'Roboto Mono', monospace; font-size: 0.68rem; color: #58a6ff; text-transform: uppercase; letter-spacing: 0.1em; margin-bottom: 12px; }
      .form-select, .form-control { background: #161b22 !important; border: 1px solid #30363d !important; color: #e6edf3 !important; border-radius: 6px !important; font-family: 'Roboto Mono', monospace !important; font-size: 0.8rem !important; }
      .form-label { font-family: 'Roboto Mono', monospace; font-size: 0.7rem; color: #8b949e; text-transform: uppercase; letter-spacing: 0.07em; }
      .btn-convert { background: linear-gradient(135deg, #1f6feb, #388bfd); border: none; border-radius: 8px; color: white; font-family: 'Roboto Mono', monospace; font-weight: 700; font-size: 0.88rem; letter-spacing: 0.05em; padding: 10px 24px; width: 100%; margin-top: 16px; transition: all 0.2s ease; }
      .btn-convert:hover { background: linear-gradient(135deg, #388bfd, #58a6ff); box-shadow: 0 0 20px #58a6ff33; transform: translateY(-1px); }
      .btn-convert:disabled { opacity: 0.4; transform: none; box-shadow: none; cursor: not-allowed; }
      .status-box { background: #0d1117; border: 1px solid #30363d; border-radius: 8px; padding: 12px 16px; font-family: 'Roboto Mono', monospace; font-size: 0.75rem; color: #8b949e; margin-top: 16px; min-height: 60px; line-height: 1.7; }
      .status-success { color: #3fb950; } .status-error { color: #f85149; } .status-running { color: #d29922; }
      .info-icon { display: inline-block; width: 14px; height: 14px; background: #30363d; border-radius: 50%; color: #8b949e; font-size: 0.6rem; font-weight: 700; text-align: center; line-height: 14px; cursor: pointer; margin-left: 5px; position: relative; vertical-align: middle; }
      .info-icon:hover .tooltip-text { visibility: visible; opacity: 1; }
      .tooltip-text { visibility: hidden; opacity: 0; background: #21262d; border: 1px solid #30363d; color: #8b949e; font-family: 'Roboto Mono', monospace; font-size: 0.68rem; border-radius: 6px; padding: 6px 10px; position: absolute; z-index: 999; bottom: 120%; left: 50%; transform: translateX(-50%); width: 200px; line-height: 1.5; transition: opacity 0.2s; pointer-events: none; }
    "))
  ),
  
  div(class = "header-bar",
      div(
        div(class = "app-title", "SPARSE CONVERTER"),
        div(class = "app-subtitle", "Convert Numeric Matrix to Sparse Matrix")
      ),
      div(class = "author-tag", "Debarshi Chakraborty, Somak Dutta, Vivekananda Roy")
  ),
  
  div(style = "padding: 0 24px 24px;",
      div(class = "step-card",
          div(class = "step-badge", "MATRIX CONVERSION"),
          div(class = "step-title", "Convert Numeric Genotype File to Sparse Matrix"),
          div(class = "step-desc",
              "Upload a numeric genotype file (should be gzipped) and convert it to a sparse matrix (.rds) saved in your chosen directory."),
          
          div(class = "upload-area",
              div(class = "upload-label", "Genotype File (.gz)"),
              fileInput("geno_file", label = NULL,
                        accept = ".gz",
                        placeholder = "genotype_file.txt.gz",
                        buttonLabel = "Browse")
          ),
          
          div(class = "param-section",
              div(class = "param-title", "⚙ CONVERSION PARAMETERS"),
              div(
                tags$label(class = "form-label", "Number of Genotypes",
                           tags$span(class = "info-icon", "i",
                                     tags$span(class = "tooltip-text", "Maximum number of genotypes (rows) to read. An upper bound is fine."))),
                numericInput("num_genotypes", label = NULL, value = 1000, min = 1, step = 1)
              ),
              div(
                tags$label(class = "form-label", "Separator",
                           tags$span(class = "info-icon", "i",
                                     tags$span(class = "tooltip-text", "Character that separates entries in each line of the file."))),
                selectInput("separator", label = NULL,
                            choices = c("Tab (\\t)" = "\t", "Comma (,)" = ",", "Space ( )" = " ", "Semicolon (;)" = ";"),
                            selected = "\t")
              ),
              div(
                tags$label(class = "form-label", "Output Filename"),
                textInput("out_filename", label = NULL, value = "sparse_matrix.rds", placeholder = "sparse_matrix.rds")
              ),
              div(style = "margin-top: 12px;",
                  tags$label(class = "form-label", "Save Location",
                             tags$span(class = "info-icon", "i",
                                       tags$span(class = "tooltip-text", "Choose where to save the .rds file. Falls back to working directory if the selected location does not exist."))),
                  selectInput("save_location", label = NULL,
                              choices = c("Downloads" = "downloads", "Desktop" = "desktop", "Working Directory" = "working"),
                              selected = "downloads"),
                  uiOutput("resolved_dir_display")
              )
          ),
          
          actionButton("convert_btn", "CONVERT TO SPARSE", class = "btn-convert", icon = icon("compress")),
          div(class = "status-box", uiOutput("convert_status"))
      )
  )
)

server <- function(input, output, session) {
  
  resolve_dir <- reactive({
    choice <- input$save_location
    home <- path.expand("~")
    candidate <- switch(choice,
                        "downloads" = file.path(home, "Downloads"),
                        "desktop"   = file.path(home, "Desktop"),
                        "working"   = getwd()
    )
    if(dir.exists(candidate)) candidate else getwd()
  })
  
  output$resolved_dir_display <- renderUI({
    tags$div(style = "font-family:'Roboto Mono',monospace; font-size:0.72rem; color:#484f58; margin-top:5px;",
             paste0("→ ", resolve_dir()))
  })
  
  observeEvent(input$convert_btn, {
    req(input$geno_file)
    
    shinyjs::disable("convert_btn")
    shinyjs::html("convert_btn", "<i class='fa fa-compress'></i> Converting...")
    on.exit({
      shinyjs::enable("convert_btn")
      shinyjs::html("convert_btn", "<i class='fa fa-compress'></i> CONVERT TO SPARSE")
    })
    
    output$convert_status <- renderUI({
      tags$span(class = "status-running", "⏳ Converting — this may take a while for large files...")
    })
    
    tryCatch({
      file_path     <- input$geno_file$datapath
      num_genotypes <- as.integer(input$num_genotypes)
      separator     <- input$separator
      out_filename  <- trimws(input$out_filename)
      if(!grepl("\\.rds$", out_filename, ignore.case = TRUE)){
        out_filename <- paste0(out_filename, ".rds")
      }
      
      X <- dense2sparse(file.name = file_path,
                        num.genotypes = num_genotypes,
                        separator = separator,
                        progress = FALSE)
      
      save_path <- file.path(resolve_dir(), out_filename)
      saveRDS(X, file = save_path)
      
      output$convert_status <- renderUI({
        tags$span(class = "status-success",
                  paste0("✔ Conversion complete!"), tags$br(),
                  paste0("Dimensions: ", nrow(X), " genotypes × ", ncol(X), " SNPs"), tags$br(),
                  paste0("Saved → ", save_path))
      })
      
    }, error = function(e) {
      output$convert_status <- renderUI({
        tags$span(class = "status-error", paste("✖ Error:", e$message))
      })
    })
  })
}

shinyApp(ui = ui, server = server)