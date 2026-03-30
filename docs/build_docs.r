# Renders vignettes and copies output to docs/ for Github Pages

rmarkdown::render(
  input = "vignettes/getting-started.Rmd",
  output_format = "rmarkdown::html_vignette",
  output_dir = "docs"
)

message("Done building vignettes.")
