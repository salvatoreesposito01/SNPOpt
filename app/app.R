# app/app.R
library(shiny)
options(shiny.maxRequestSize = 1024 * 1024^2) # 1 GB upload limit
source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)
