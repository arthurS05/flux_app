# app.R
# Shiny application for processing LI-COR flux data (CH₄ and CO₂)
# -------------------------------------------------------------------
# Instructions:
# - Place this file at the root of your GitHub repository.
# - Put raw CSV files (all_data.csv and Source05.csv) into a folder named "data" in the project root.
# - Install required packages: shiny, ggplot2, DT, here

# 1) Load libraries ----
library(shiny)
library(ggplot2)
library(DT)
library(here)

# 2) Define data folder ----
data_dir <- here::here("data")

# 3) File paths ----()
all_data_path <- file.path(data_dir, "all_data.csv")
list_path     <- file.path(data_dir, "list_all.csv")

# 4) Read raw CSVs ----
df_data_init <- read.csv(all_data_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
list_df0     <- read.csv(list_path,     header = TRUE, sep = ";", stringsAsFactors = FALSE)

list_df0 <- list_df0[list_df0$Type == "Flux", ]
# 5) Subset/clean your list if needed ----
list_df <- list_df0
# e.g. list_df <- list_df0[, -c(9:12)]

# 6) Experiment constants (USER MUST SET) ----
# These depend on your specific chamber setup:
#   • voltotal: total gas volume inside chamber (in liters)
#   • surface_area: chamber’s measured surface area (in m²)
voltotal     <- 11.60365     # user-defined, e.g. 11.60 L
surface_area <- 0.085633564  # user-defined, e.g. π*(0.12 m)² ≈ 0.0856 m²

#--------------------------------------------------------------------------------
# UI definition
#--------------------------------------------------------------------------------
ui <- fluidPage(
  # Application title
  titlePanel(
    strong("LI-COR Flux Data Processing", style = "font-variant:small-caps; color:#FFFFFF; background-color:#16697A; padding:10px;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      # Button to advance through the metadata list
      actionButton("nextButton", "Next sample", icon = icon("arrow-right")),
      hr(),
      # Settings for CH4 window selection
      h4("CH₄: Seconds before and after reference point"),
      numericInput("before_CH4", "Seconds before (CH₄)", value = 0),
      numericInput("after_CH4",  "Seconds after (CH₄)",  value = 300),
      numericInput("inter_CH4",  "Zoom-out padding (CH₄)", value = 20, min = 0),
      hr(),
      # Settings for CO2 window selection
      h4("CO₂: Seconds before and after reference point"),
      numericInput("before_CO2", "Seconds before (CO₂)", value = 0),
      numericInput("after_CO2",  "Seconds after (CO₂)",  value = 300),
      numericInput("inter_CO2",  "Zoom-out padding (CO₂)", value = 20, min = 0),
      hr(),
      # Button to save computed fluxes
      actionButton("button_save", "Save computed fluxes", icon = icon("floppy-disk"))
    ),
    
    mainPanel(
      # Display current metadata row and summary
      verbatimTextOutput("line"),
      verbatimTextOutput("line2"),
      
      # Tabs for CH4 and CO2 analysis
      tabsetPanel(
        tabPanel("CH₄",
                 h3("CH₄ Time Series"),
                 plotOutput("plot_CH4", click = "plot_click_CH4"),
                 verbatimTextOutput("info1_CH4"),  # Regression summary
                 verbatimTextOutput("info2_CH4")   # Flux calculation summary
        ),
        tabPanel("CO₂",
                 h3("CO₂ Time Series"),
                 plotOutput("plot_CO2", click = "plot_click_CO2"),
                 verbatimTextOutput("info1_CO2"),
                 verbatimTextOutput("info2_CO2")
        )
      ),
      
      hr(),
      # Final results table
      DT::dataTableOutput("data_tablefin")
    )
  )
)

#--------------------------------------------------------------------------------
# Server logic
#--------------------------------------------------------------------------------
server <- function(input, output, session) {
  # Reactive counter to track current row in list_df
  counter <- reactiveVal(0)
  maxLine <- nrow(list_df)
  
  # Prepare the main data frame: drop header lines and rename columns
  df0 <- df_data_init[-c(1:4,6), ]              # remove metadata rows
  colnames(df0) <- df0[1, ]                     # set column names
  df_main <- df0[-1, ]                          # drop the row used for names
  
  # Extract CH4 and CO2 time series into separate data frames
  df_CH4 <- df_main[, c("DATE", "TIME", "CH4")]
  df_CH4$CH4 <- as.numeric(df_CH4$CH4)
  
  df_CO2 <- df_main[, c("DATE", "TIME", "CO2")]
  df_CO2$CO2 <- as.numeric(df_CO2$CO2)
  
  # Reactive storage for computed results
  results <- reactiveValues(
    data = data.frame(
      Files = character(), Lake = character(), Type = character(),
      Name_sample = character(), remark = character(), Name_LICOR = character(),
      date_list = character(), H_list = character(), Gaz = character(),
      magnitude_value = character(), H_begin = character(), H_end = character(),
      slope = numeric(), Rsquare = numeric(), volume = numeric(),
      surface = numeric(), Flux_mg_m2_s = numeric(),
      stringsAsFactors = FALSE
    )
  )
  
  #––– NEW: hold the last regression/flux results for each gas –––––––––––––––––
  last_CH4 <- reactiveValues(
    begin = NULL, end = NULL, slope = NULL, Rsq = NULL, flux = NULL
  )
  last_CO2 <- reactiveValues(
    begin = NULL, end = NULL, slope = NULL, Rsq = NULL, flux = NULL
  )
  
  # Helper function to compute slope and R² given a vector of concentrations
  compute_regression <- function(values) {
    n <- length(values)
    df <- data.frame(time_index = seq_len(n), conc = values)
    fit <- lm(conc ~ time_index, data = df)
    corr <- cor.test(df$conc, df$time_index)
    list(
      slope = coef(fit)[2],
      Rsquare = corr$estimate^2
    )
  }
  
  # Observe "Next sample" button to advance through metadata list
  observeEvent(input$nextButton, {
    # Increment counter (but do not exceed maxLine)
    new_idx <- min(counter() + 1, maxLine)
    counter(new_idx)
    
    # Extract metadata for current sample
    row <- list_df[new_idx, ]
    
    # Display raw metadata
    output$line <- renderPrint({ row })
    output$line2 <- renderPrint({
      paste0("Sample ", new_idx, "/", maxLine,
             " : date=", row$date,
             " time=", row$heure_debut)
    })
    
    ## CH₄: plot and regression ##
    output$plot_CH4 <- renderPlot({
      idx0 <- which(df_CH4$DATE == row$date & df_CH4$TIME == row$heure_debut)
      before <- input$before_CH4; after <- input$after_CH4; pad <- input$inter_CH4
      low_reg  <- idx0 - before
      high_reg <- idx0 + after - 1
      low_plot <- low_reg - pad
      high_plot<- high_reg + pad
      
      plot_data <- df_CH4[low_plot:high_plot, , drop = FALSE]
      plot_data$phase <- rep(
        c("padding", "analysis", "padding"),
        times = c(pad, nrow(plot_data) - 2*pad, pad)
      )
      
      ggplot(plot_data, aes(x = TIME, y = CH4, color = phase)) +
        geom_point() +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "CH₄ Flux Regression",
             x = "Time", y = "Concentration (ppb)") +
        theme_light()
    })
    
    output$info1_CH4 <- renderPrint({
      idx0 <- which(df_CH4$DATE == row$date & df_CH4$TIME == row$heure_debut)
      before <- input$before_CH4; after <- input$after_CH4
      reg_vals <- df_CH4[(idx0 - before):(idx0 + after - 1), "CH4"]
      res <- compute_regression(reg_vals)
      cat(
        "Start time:", df_CH4$TIME[idx0 - before], "\n",
        "End time:  ", df_CH4$TIME[idx0 + after - 1], "\n",
        sprintf("Slope: %.4f ppb/s\nR²: %.4f", res$slope, res$Rsquare)
      )
    })
    
    output$info2_CH4 <- renderPrint({
      idx0 <- which(df_CH4$DATE == row$date & df_CH4$TIME == row$heure_debut)
      before <- input$before_CH4; after <- input$after_CH4
      reg_vals <- df_CH4[(idx0 - before):(idx0 + after - 1), "CH4"]
      res <- compute_regression(reg_vals)
      
      # Convert slope (ppb/s) → flux (mg/m2/s)
      slope_mol_L_s <- (res$slope * 1e-9) / 23
      slope_mol_s   <- slope_mol_L_s * voltotal
      flux_mol_m2_s <- slope_mol_s / surface_area
      flux_mg_m2_s  <- flux_mol_m2_s * 16.04 * 1000
      
      # Store for saving
      last_CH4$begin <- df_CH4$TIME[idx0 - before]
      last_CH4$end   <- df_CH4$TIME[idx0 + after - 1]
      last_CH4$slope <- res$slope
      last_CH4$Rsq   <- res$Rsquare
      last_CH4$flux  <- flux_mg_m2_s
      
      cat(sprintf("Computed CH₄ flux: %.3f mg/m2/s", flux_mg_m2_s))
    })
    
    ## CO₂: same workflow as CH₄ ##
    output$plot_CO2 <- renderPlot({
      idx0 <- which(df_CO2$DATE == row$date & df_CO2$TIME == row$heure_debut)
      before <- input$before_CO2; after <- input$after_CO2; pad <- input$inter_CO2
      low_reg  <- idx0 - before
      high_reg <- idx0 + after - 1
      low_plot <- low_reg - pad
      high_plot<- high_reg + pad
      
      plot_data <- df_CO2[low_plot:high_plot, , drop = FALSE]
      plot_data$phase <- rep(
        c("padding", "analysis", "padding"),
        times = c(pad, nrow(plot_data) - 2*pad, pad)
      )
      
      ggplot(plot_data, aes(x = TIME, y = CO2, color = phase)) +
        geom_point() +
        scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "CO₂ Flux Regression",
             x = "Time", y = "Concentration (ppm)") +
        theme_light()
    })
    
    output$info1_CO2 <- renderPrint({
      idx0 <- which(df_CO2$DATE == row$date & df_CO2$TIME == row$heure_debut)
      before <- input$before_CO2; after <- input$after_CO2
      reg_vals <- df_CO2[(idx0 - before):(idx0 + after - 1), "CO2"]
      res <- compute_regression(reg_vals)
      cat(
        "Start time:", df_CO2$TIME[idx0 - before], "\n",
        "End time:  ", df_CO2$TIME[idx0 + after - 1], "\n",
        sprintf("Slope: %.4f ppm/s\nR²: %.4f", res$slope, res$Rsquare)
      )
    })
    
    output$info2_CO2 <- renderPrint({
      idx0 <- which(df_CO2$DATE == row$date & df_CO2$TIME == row$heure_debut)
      before <- input$before_CO2; after <- input$after_CO2
      reg_vals <- df_CO2[(idx0 - before):(idx0 + after - 1), "CO2"]
      res <- compute_regression(reg_vals)
      
      # Convert slope (ppm/s) → flux (mg/m2/s)
      slope_mol_L_s <- (res$slope * 1e-6) / 23
      slope_mol_s   <- slope_mol_L_s * voltotal
      flux_mol_m2_s <- slope_mol_s / surface_area
      flux_mg_m2_s  <- flux_mol_m2_s * 44.01 * 1000
      
      # Store for saving
      last_CO2$begin <- df_CO2$TIME[idx0 - before]
      last_CO2$end   <- df_CO2$TIME[idx0 + after - 1]
      last_CO2$slope <- res$slope
      last_CO2$Rsq   <- res$Rsquare
      last_CO2$flux  <- flux_mg_m2_s
      
      cat(sprintf("Computed CO₂ flux: %.3f mg/m2/s", flux_mg_m2_s))
    })
  })
  
  # Observe Save button to store current CH4 and CO2 results
  observeEvent(input$button_save, {
    req(counter() > 0)  # ensure a sample is selected
    
    # Base metadata
    row   <- list_df[counter(), ]
    base  <- list(
      Files       = row$Fichier_de_donnees,
      Lake        = row$lake,
      Type        = row$Type,
      Name_sample = row$Nom_echantillon_BD,
      remark      = row$remarque,
      Name_LICOR  = row$Nom_dans_fichier_LICOR,
      date_list   = row$date,
      H_list      = row$heure_debut
    )
    
    # CH4 entry
    ch4_row <- data.frame(
      base,
      Gaz             = "CH4",
      magnitude_value = "ppb",
      H_begin         = last_CH4$begin,
      H_end           = last_CH4$end,
      slope           = last_CH4$slope,
      Rsquare         = last_CH4$Rsq,
      volume          = voltotal,
      surface         = surface_area,
      Flux_mg_m2_s    = last_CH4$flux,
      stringsAsFactors = FALSE
    )
    # CO2 entry
    co2_row <- data.frame(
      base,
      Gaz             = "CO2",
      magnitude_value = "ppm",
      H_begin         = last_CO2$begin,
      H_end           = last_CO2$end,
      slope           = last_CO2$slope,
      Rsquare         = last_CO2$Rsq,
      volume          = voltotal,
      surface         = surface_area,
      Flux_mg_m2_s    = last_CO2$flux,
      stringsAsFactors = FALSE
    )
    
    # Append rows
    results$data <- rbind(results$data, ch4_row, co2_row)
    assign("df_flux_final", isolate(results$data), envir = .GlobalEnv)
    
    # Re-render table
    output$data_tablefin <- DT::renderDataTable({
      DT::datatable(results$data, options = list(pageLength = 10))
    })
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
