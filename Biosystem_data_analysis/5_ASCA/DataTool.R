library(shiny)
library(shinydashboard)
library(tidyverse)
library(plyr)
library(ggrepel)


source("DataInput.R")
source("Preprocess.R")
source("ASCA_module.R")
source("UniVariate.R")

ui <- dashboardPage(
  dashboardHeader(title = "Data analysis tool",
    dropdownMenu(type = "notifications",
      notificationItem(text = "Warning",
                       icon = icon("exclamation-triangle"),
                       status = "warning")
                )
  ),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Data input", tabName = "datainput", icon = icon("file")),
      menuItem("Preprocess", tabName = "preprocess", icon = icon("steam-square")),
      menuItem("Univariate", tabName = "univariate", icon = icon("steam-square")),
      menuItem("ASCA", tabName = "asca", icon = icon("steam-square"),
        menuSubItem("Build Model", tabName = "buildmodel", icon = icon("steam-square")),
        
        menuItem("Plots", tabName = "plots", icon = icon("image"),
          menuSubItem("Factors", tabName = "factors", icon = icon("image")),
          menuSubItem("Interactions", tabName = "interactions", icon = icon("image")),
          menuSubItem("Combination", tabName = "combinations", icon = icon("image")),
          menuSubItem("Residuals", tabName = "residuals", icon = icon("image"))
        )
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem("datainput", readFilesUI(id = "File")
      ),
      tabItem("preprocess", preProcessUI(id = "Preprocess")
      ),
      tabItem("univariate", uniVariateUI(id = "Univariate")
      ),
      tabItem("buildmodel", buildModelUI(id = "ASCA")
      ),
      tabItem("factors", makeFactorsUI(id = "ASCA")  
      ),
      tabItem("interactions", makeInteractionsUI(id = "ASCA")  
      ), 
      tabItem("combinations", makeCombinationsUI(id = "ASCA") 
      ),
      tabItem("residuals", makeResidualsUI(id = "ASCA") 
      )
    )
  )
)

server <- function(input, output, session) {
  # DAT.session is used for variables needed by other modules
  DAT.session <- reactiveValues()
  
  # the variable data is used for some variables instead of DAT.session because
  # they are only used in preProcess.
  data <- callModule(readFiles, id = "File", DAT.session)
  callModule(preProcess, id = "Preprocess", data, DAT.session)
  callModule(uniVariate, id = "Univariate", DAT.session)
  callModule(ASCAmodel, id = "ASCA", DAT.session)
}

shinyApp(ui, server)