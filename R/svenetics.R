#' @title Launch the SVENETICS Shiny App
#' @description Opens the SVENETICS GUI in your browser.
#' @export
svenetics <- function() {
  app_dir <- system.file("app", package = "bravo")
  if (app_dir == "") {
    stop("Could not find the app directory. Please try re-installing svenetics.")
  }
  runApp(app_dir, display.mode = "normal")
}

