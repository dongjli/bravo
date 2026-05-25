#' @title Launch the Sparse Converter Shiny App
#' @description Opens the Sparse Matrix Converter GUI in your browser.
#' @export
dense_to_sparse_converter <- function() {
  app_dir <- system.file("app2", package = "bravo")
  if (app_dir == "") {
    stop("Could not find the app directory. Try re-installing bravo.")
  }
  runApp(app_dir, display.mode = "normal")
}