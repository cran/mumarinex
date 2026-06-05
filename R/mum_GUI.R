#' mum_GUI
#'
#' @title Online Graphical User Interface to compute MUMARINEX
#'
#' @description
#' Launch the online Graphical User Interface to compute the MUMARINEX index and its sub-indices.
#'
#' @returns Opens the web browser to the online Shiny application.
#' @export
#'
#' @importFrom utils browseURL
#'
#' @examples
#' mum_GUI()
mum_GUI<-function() {
  browseURL("https://mumarinex.shinyapps.io/mumarinex_online_interface/")
  message("Opening the MUMARINEX GUI in your browser")
}
