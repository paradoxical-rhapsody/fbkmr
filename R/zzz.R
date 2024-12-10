.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Please only use skmbayes function for now.")
  packageStartupMessage("The package requires Gurobi package. Please check the website for installation: https://www.gurobi.com/.")
}

release_questions <- function() {
  c(
    "Have you updated the vignette and posted to GitHub?"
  )
}
