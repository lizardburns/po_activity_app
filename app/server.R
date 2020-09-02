########### PO activity calculator - v1.2 - server ============================

server <- function(input, output) {
  # Prepare absorbance results file -------------------------------------------
  getAbsorbance <- reactive({
    out <- list()
    
    if(is.null(input$absorbResults)) {
      out <- NULL
    } else if(input$wavelength == "All three (475nm, 600nm, 490nm)" & !is.null(input$absorbResults)){
      inFile <- input$absorbResults
      dabsorb <- read.xlsx(inFile$datapath,
                           sheet = 1, 
                           colNames = FALSE, 
                           rowNames = FALSE)
      
      starts <- str_which(str_to_lower(dabsorb$X1), "cycle")
      data_names <- dabsorb[starts[1],]
      
      out$l475 <- dabsorb %>% 
        slice(starts[1]+1:13) %>% 
        mutate(across(everything(), as.numeric))
      out$l600 <- dabsorb %>% 
        slice(starts[2]+1:13) %>% 
        mutate(across(everything(), as.numeric))
      out$l490 <- dabsorb %>% 
        slice(starts[3]+1:13) %>% 
        mutate(across(everything(), as.numeric))
      
      names(out$l475) <- data_names
      names(out$l600) <- data_names
      names(out$l490) <- data_names
      
     } else if(input$wavelength == "490nm" & !is.null(input$absorbResults)){
      # this table need transposed
      inFile <- input$absorbResults
      dabsorb <- read.xlsx(inFile$datapath,
                           sheet = 1,
                           colNames = FALSE, rowNames = FALSE)
      
      starts <- str_which(str_to_lower(dabsorb$X1), "cycle")
      l490 <- dabsorb %>% slice(starts:nrow(dabsorb))
      rownames(l490) <- l490$X1
      l490 <- l490 %>% mutate(across(-X1, as.character))
      
      out$l490 <- l490 %>% 
        tidyr::pivot_longer(-X1) %>% 
        tidyr::pivot_wider(names_from = X1, values_from = value) %>% 
        dplyr::select(-name) %>% 
        mutate(across(everything(), as.numeric))
      
     }
    return(out)
  })
  
  # Prepare protein results file ----------------------------------------------
  getProtein <- reactive({
    wells <- c("A1", "A3", "A5", "A7", "A9", "A11", "B1", "B3", "B5", "B7", "B9", "B11", "C1", "C3", "C5", "C7", "C9", "C11", "D1", "D3", "D5", "D7", "D9", "D11", "E1", "E3", "E5", "E7", "E9", "E11", "F1", "F3", "F5", "F7", "F9", "F11", "G1", "G3", "G5", "G7", "G9", "G11", "H1", "H3", "H5", "H7", "H9", "H11")
    
    if(is.null(input$proteinResults)) return(NULL)
    else {
      proteinFile <- input$proteinResults
      dprotein <- read_csv(proteinFile$datapath) %>% 
        mutate(well = str_to_upper(well)) %>% 
        filter(well %in% wells)
      
      if(str_detect(str_to_lower(paste(dprotein$sampleID, collapse = " ")), "blank")){
        dprotein$mg_protein[str_to_lower(dprotein$sampleID) == "blank"] <- 1
      }
      
      return(dprotein)
    }
  })
  
  # Signal to ready UI for results download -----------------------------------
  output$filesUploaded <- reactive({
    if(is.null(getAbsorbance()) | is.null(getProtein())){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  outputOptions(output, 'filesUploaded', suspendWhenHidden=FALSE)
  
  # Blank needed? -------------------------------------------------------------
  # get user input if case insensitive string "blank" can't be found in the 
  # protein file containing sample IDs
  output$blankneeded <- reactive({
    if(is.null(getAbsorbance()) | is.null(getProtein())){
      return(FALSE)
    } else if(!is.null(getAbsorbance()) & !is.null(getProtein()) & str_detect(str_to_lower(paste(getProtein()$sampleID, collapse = " ")), "blank")){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  outputOptions(output, 'blankneeded', suspendWhenHidden=FALSE)
  
  # Analysis ------------------------------------------------------------------
  plate_calcs475 <- reactive({
    validate(need(!is.null(getAbsorbance()), "Waiting for file upload"))
    if(input$wavelength == "All three (475nm, 600nm, 490nm)"){
      calculate_dabs_plots(getAbsorbance()$l475)
    } else {
      NULL
    }
  })
  
  plate_calcs490 <- reactive({
    validate(need(!is.null(getAbsorbance()), "Waiting for file upload"))
    calculate_dabs_plots(getAbsorbance()$l490)
  })
  
  plate_calcs600 <- reactive({
    validate(need(!is.null(getAbsorbance()), "Waiting for file upload"))
    if(input$wavelength == "All three (475nm, 600nm, 490nm)"){
      calculate_dabs_plots(getAbsorbance()$l600)
    } else {
      NULL
    }
  })
  
  dabs475 <- reactive({
    validate(need(!is.null(plate_calcs475()), "No velocities available for lambda 475"))
    map_dbl(plate_calcs475(), ~.x$vmax)
    })
  
  dabs600 <- reactive({
    validate(need(!is.null(plate_calcs600()), "No velocities available for lambda 600"))
    map_dbl(plate_calcs600(), ~.x$vmax)
  })
  
  dabs490 <- reactive({
    validate(need(!is.null(plate_calcs490()), "No velocities available for lambda 490"))
    map_dbl(plate_calcs490(), ~.x$vmax)
  })
  
  plots475 <- reactive({
    # validate(need(!is.null(plate_calcs475()), "No velocities available for lambda 475"))
    if(is.null(plate_calcs475())){
      NULL
    } else {
      map(plate_calcs475(), ~.x$p)
    }
  })
  
  plots600 <- reactive({
    # validate(need(!is.null(plate_calcs600()), "No velocities available for lambda 600"))
    if(is.null(plate_calcs600())){
      NULL
    } else {
      map(plate_calcs600(), ~.x$p)
    }
  })
  
  plots490 <- reactive({
    # validate(need(!is.null(plate_calcs490()), "No velocities available for lambda 490"))
    if(is.null(plate_calcs490())){
      NULL
    } else {
      map(plate_calcs490(), ~.x$p)
    }
  })
  
  # generate results ----------------------------------------------------------
  results <- reactive({
    bval475 <- input$blank475
    bval490 <- input$blank490
    
    if(!is.null(getAbsorbance()) & !is.null(getProtein())){
      if(input$wavelength == "All three (475nm, 600nm, 490nm)"){
        calculate_results_summary(lambda475 = dabs475(), lambda600 = dabs600(), lambda490 = dabs490(), 
                                  protein_table = getProtein(), 
                                  blank475 = bval475, blank490 = bval490)
      } else if(input$wavelength == "490nm"){
        calculate_490_summary(dabs490(), getProtein(), blank490 = bval490)
      }
    } else {
      NULL
    }
  })
  
  # prepare outputs -----------------------------------------------------------
  observe({
    if(!is.null(results())) write_csv(results(), "results.csv")
  })
  
  output$results_dload <- downloadHandler(
    filename = paste0(now(), "_results.csv"),
    content = function(file){file.copy("results.csv", file)}
  )
  
  output$plot1 <- renderPlotly({
    validate(
      need(!is.null(results()), "Results will be plotted here once data is uploaded")
    )
    
    if(input$wavelength == "All three (475nm, 600nm, 490nm)"){
      p <- ggplot(results(), aes(x = sampleID, y = mean_dabs_min_475_600_mgprotein, colour = mean_dabs_min_475_600_mgprotein)) +
        geom_point()
    } else if(input$wavelength == "490nm"){
      p <- ggplot(results(), aes(x = sampleID, y = mean_dabs_min_490_mgprotein, colour = mean_dabs_min_490_mgprotein)) +
        geom_point()
    }
    
    p <- p +
      scale_colour_viridis() +
      theme_bw() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
      labs(y = "PO activity /mg/min",
           colour = "PO activity /mg/min",
           title = "Results summary")
    ggplotly(p)
  })
  
  output$report_dload <- downloadHandler(
    filename = "report.html",
    content = function(file) {
      withProgress(message = 'Preparing report, please wait!', {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list("plots_lambda475" = plots475(),
                       "plots_lambda600" = plots600(),
                       "plots_lambda490" = plots490(),
                       rendered_by_shiny = TRUE)

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in the app).
        rmarkdown::render(
          input = tempReport, 
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
      })
    }
  )
  
}
