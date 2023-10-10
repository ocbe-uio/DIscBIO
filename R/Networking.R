#' @title Plotting the network.
#' @description This function uses STRING API to plot the network.
#' @references https://string-db.org/api/
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the
#' resulted network. If `NULL`, the network will be saved to a temporary
#' directory
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @param retries maximum number of attempts to connect to the STRING api.
#' @importFrom utils download.file
#' @importFrom png readPNG
#' @importFrom graphics plot rasterImage
#' @return A plot of the network
Networking <- function(
    data,
    FileName = NULL,
    species = "9606",
    plot_width = 25,
    plot_height = 15,
    retries = 3) {
  # ======================================================== #
  # Validation                                               #
  # ======================================================== #
  if (length(data) > 600) stop("Your gene list is too big")

  # ======================================================== #
  # Processing                                               #
  # ======================================================== #
  repos <- retrieveURL(data, species, "highres_image")
  y <- repos$request$url
  if (!is.null(FileName)) {
    FileName <- paste0("network", FileName, ".png")
  } else {
    FileName <- tempfile()
  }
  download.file(y, FileName, mode = "wb")
  message(
    "\n",
    "You can see the network with high resolution ",
    "by clicking on the following link:",
    "\n",
    paste0(y)
  )
  Network <- readPNG(FileName, native = TRUE)
  set_plot_dimensions <- function(width_choice, height_choice) {
    opar <- withr::with_options(
      repr.plot.width = width_choice,
      repr.plot.height = height_choice
    )
    on.exit(withr::with_options(opar))
  }
  set_plot_dimensions(plot_width, plot_height)

  # Plotting ----------------------------------------------- #
  plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
  rasterImage(Network, 0, 0, 1, 1)
  set_plot_dimensions(8, 8) # resets to default values
}
