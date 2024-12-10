.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Please only use skmbayes function for now.")
  packageStartupMessage("The package requires Gurobi package. Please check the website for installation: https://github.com/junwei-lu/fbkmr.")
}

release_questions <- function() {
  c(
    "Have you checked the Github repositorary?"
  )
}
