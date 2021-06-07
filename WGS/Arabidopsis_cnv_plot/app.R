#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

# load data
h60 <- "h60d_rmdup_Atchr_only.hits-vs-At.mapped.renamechr.hits.log2-0.6.pvalue-0.001.minw-4.cnv"
h300 <- "h300d_rmdup_Atchr_only.hits-vs-At.mapped.renamechr.hits.log2-0.6.pvalue-0.001.minw-4.cnv"
# data <- read.table("h60d_rmdup_Atchr_only.hits-vs-At.mapped.renamechr.hits.log2-0.6.pvalue-0.001.minw-4.cnv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
# h300d <- read.table("h300d_rmdup_Atchr_only.hits-vs-At.mapped.renamechr.hits.log2-0.6.pvalue-0.001.minw-4.cnv")

### possibility of doing a zoomable plot, and/or plotting h60 and h300 in top bottom fashion for comparison
# Define UI for application that draws a scatterplot
ui <- fluidPage(
  
  # Application title
  titlePanel("CNV plot"),
  
  # Sidebar with a slider input for range of chromosome position 
  sidebarLayout(
    sidebarPanel(
      sliderInput("range", "Chromosome range",
                  min = 1, max = 30000000,
                  value = c(1,1000000)),
      ## slider range
      # sliderInput("slider1", label = h3("Slider"), min = 0, 
      #               max = 100, value = 50)
      # ),
      
      ## checkbox group
      checkboxGroupInput("chromosome", label = h3("Chromosome"),
                         choices = list("Chromosome 1" = 1, "Chromosome 2" = 2, "Chromosome 3" = 3, "Chromosome 4" = 4, "Chromosome 5" = 5),
                         selected = 1),
      # select box 
      # selectInput("select", label = h3("Select box"), 
      #             choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), 
      #             selected = 1),    
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput(outputId = "cnvPlot", height = 300),
      plotOutput(outputId = "cnvPlot2", height = 300)
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  ## ggplot for cnv
  plot_cnv <- function(cnv_data) {
    cnv <- read.table(cnv_data, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    cnv <- cnv[is.finite(cnv$log2),]
    cnv <- cnv[cnv$chromosome == input$chromosome, ]
    g <- ggplot(cnv, aes(x=start, y=log2))
    g + geom_point(size=0.1) + facet_grid(chromosome ~ . ) + geom_hline(yintercept = 0, linetype = "dashed", color = "slategrey")  + labs(x = "position", y = "log2 ratio") + theme(axis.ticks = element_blank()) + xlim(input$range[1], input$range[2]) + ylim(-10, 7.5) + ggtitle(as.character(cnv_data))
  }
  output$cnvPlot <- renderPlot({
    # either try how to make facet grid work, or create another selection bar for chromosome
    # setting chromosome to be drawn
    # data <- data[is.finite(data$log2),]
    # data <- data[data$chromosome == input$chromosome, ]
    
    # plot scatterplot of cnv along chromosome
    # g <- ggplot(data,aes(x=start, y=log2))
    # g + geom_point(size=0.1)  + facet_grid(chromosome ~ . ) + geom_hline(yintercept = 0, linetype = "dashed", color = "slategrey")  + labs(x = "position", y = "log2 ratio") + theme(axis.ticks = element_blank()) + 
    #   xlim(input$range[1], input$range[2]) 
    
    
    # try calling function instead
    plot_cnv(h60)
  })
  
  output$cnvPlot2 <- renderPlot({
    plot_cnv(h300)
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

