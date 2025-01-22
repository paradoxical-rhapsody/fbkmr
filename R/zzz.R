.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Please only use skmbayes function for now.")
  packageStartupMessage("The package requires Gurobi package. Please check the website for installation: https://www.gurobi.com/.")
  packageStartupMessage("Use following function to install the package: install.packages(path-to-file, repos = NULL, type=\"source\",dependencies = TRUE, INSTALL_opts=\"--no-multiarch\")")
}

release_questions <- function() {
  c(
    "Have you updated the vignette and posted to GitHub?"
  )
}
