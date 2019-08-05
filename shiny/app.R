library(shiny)
library(ggplot2)
load('spls_out.RData')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Pairwise correlations for sPLS signature genes"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "anno", label = "Genomic context whose methylation link with Expression is studied", 
                        choices = c("Promter"="prom", "Genebody"="genebody", "Promoter-CGI"='prom_cgi', "Promoter-nonCGI"='prom_noncgi')),
            selectInput(inputId = "comp", label = "sPLS component", choices = 1:3),
            selectInput(inputId = "omic", label = "The 'omic from which the signature is shown", choices = c(Methylation="BSseq", Expression="RNAseq"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("boxplots"),
           verbatimTextOutput("text"),
           tableOutput("df")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$boxplots <- renderPlot({
        spls_out[["prom"]][1][[input$omic]][["plots"]]$X
        spls_out[[input$anno]][[paste0('comp_', input$comp)]][[input$omic]][["plots"]]
        # hist(rnorm(100), main=paste(class(input$anno), class(input$comp), input$omic, collapse = " "))
    })
      
    output$text <- renderText({
      "\n The signature table (Those with missing gene IDs came from non canonical chropmosomes which can be retrieved if needed): \n"
    })
    
    output$df <- renderTable({
        spls_out[[input$anno]][[paste0('comp_', input$comp)]][[input$omic]][["signature"]]
        # hist(rnorm(100), main=paste(class(input$anno), class(input$comp), input$omic, collapse = " "))
    }, rownames = TRUE)  
}

# Run the application 
shinyApp(ui = ui, server = server)
