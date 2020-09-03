########### PO activity calculator - v1.2 - ui ================================

ui <- fluidPage(
  column(2),
  column(8,
         titlePanel("PO Activity calculator"),
         br(),
         helpText("Given absorbance assay data from a 96 well plate (in xlsx format) and corresponding protein mass data and sample names linked to absorbance data by well id (in csv format), will calculate max velocity and enzyme activity standardised by protein mass."),
         br(),
         radioButtons("wavelength", "Which wavelengths do you have data for?", choices = c("490nm", "All three (475nm, 600nm, 490nm)")),
         br(),
         helpText(strong("Analysis and preparation of results file will commence once files are uploaded.")),
         br(),
         helpText("Browse your computer for absorbance assay results (xlsx file)"),
         fileInput("absorbResults", "Navigate to xlsx file:"),
         br(),
         helpText("Browse your computer for protein mass data (csv file)"),
         fileInput("proteinResults", "Navigate to csv file:"),
         br(),
         conditionalPanel(
           condition = "output.blankneeded",
           helpText(strong("There does not appear to be a blank on the plate.")),
           helpText("Please supply the absorbance value to adjust mean absorbance of each sample by."),
           helpText("Type a numeric value for the relevant blank."),
           numericInput("blank490", "490nm blank: input a value", min = -1, max = 1, value = 0), 
           conditionalPanel(
             condition = 'input.wavelength == "All three"',
             numericInput("blank475", "475-600nm blank: input a value", min = -1, max = 1, value = 0)
           )
         ),
         br(),
         conditionalPanel(
           condition = "output.filesUploaded",
           helpText("Your results and report (with absorbance vs time plots) are available to download by clicking the buttons below. Results are also plotted below."),
           br(),
           fluidRow(column(3, downloadButton("results_dload", "Download results")),
                    column(2),
                    column(3, downloadButton("report_dload", "Download report"))),
           br(),
           helpText("Refresh the browser to analyse another assay run or view results summary below."),
           br(),
           helpText("Hover over points to reveal the sample ID"),
           plotlyOutput("plot1")
         ),
  ),
  column(2)
)
