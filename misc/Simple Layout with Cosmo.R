library(shinythemes)

server <- function(input, output, session) {
  
}

ui <- fluidPage(theme=shinytheme("cosmo"),
                
                titlePanel("Shiny Theme is cosmo"),
                
                sidebarLayout(
                  
                  sidebarPanel(
                    h3("Upload countMatrix"),
                    actionButton("button", "Upload")
                  ), 
                  
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Heatmap"), 
                      tabPanel("Volcano Plot"), 
                      tabPanel("Count Matrix")
                    )
                  )
                )
)

shinyApp(ui = ui, server = server)
