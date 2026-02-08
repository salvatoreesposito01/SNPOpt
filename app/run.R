# run.R — launch SNPoptimizer Shiny app
library(shiny)
options(shiny.maxRequestSize = 1024 * 1024^2) # 1 GB

app_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

# Se hai ui.R/server.R oppure app.R nella root, va bene.
# Se invece li tieni in una sottocartella "app", cambia a "app".
runApp(app_dir, launch.browser = TRUE)
