#!/usr/bin/Rscript
#  GPArcLen/app.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.08.2018

library(shiny)


set.seed(123)
n <- 10
X <- matrix(rnorm(2*n), ncol = 2)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
    includeScript("www/fun.js"),
    # App title ----
    titlePanel("Evaluating GP Stationarity"),
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        # make_inputs creates a sideBarPanel object with the appropriate inputs
        actionButton("sis", "Split Input Space"),
        # Main panel for displaying outputs ----
         mainPanel(scatterD3Output("distPlot"))
    )
)

# Server logic to fit and visualize the model with the user specified hyperparams
server <- function(input, output) {
  output$distPlot <- renderScatterD3({
    #A JS callback func
    lasso_callback <- "
    function(sel) {
        //alert(sel.data().map(function(d) {return d.lab}).join('\\n'));
        d = sel.data();
    }
    "

        print(input$meme_dankness)
        scatterD3(x = X[,1], y = X[,2], lasso = TRUE, lasso_callback = lasso_callback)
  })
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
