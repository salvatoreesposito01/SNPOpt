# install.R — install dependencies for SNPoptimizer
repos <- "https://cloud.r-project.org"

pkgs <- c(
  "shiny",
  "GA",
  "data.table",
  "ggplot2"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) {
  install.packages(to_install, repos = repos)
}

cat("✅ Done. Installed packages:\n")
print(pkgs)
