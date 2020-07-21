#' @title Plotting the network.
#' @description This function uses STRING-api to plot the network.
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the
#' resulted network. If `NULL`, the network will be saved to a temporary
#' directory
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @importFrom httr GET status_code
#' @importFrom utils download.file
#' @importFrom png readPNG
#' @importFrom graphics plot rasterImage
#' @return A plot of the network
Networking <- function(
    data,
    FileName    = NULL,
    species     = "9606",
    plot_width  = 25,
    plot_height = 15
    )
    {
        if (length(data) > 600) {
            print("Your gene list is too big")
        } else{
            string_api_url <- "https://string-db.org/api/"
            output_format <- "highres_image"
            method <- "network"
            your_identifiers <- ""
            optional_parameters <- ""

            # Construct API request
            genes <- data
            repos <-
                GET(
                    url = paste0(
                        string_api_url,
                        output_format,
                        '/',
                        method,
                        '?identifiers=',
                        paste(as.character(data), collapse = "%0d"),
                        "&",
                        "species=",
                        species
                    )
                )
            message(
                "Examine response components =",
                status_code(repos),
                "\t",
                "(200 means successful)",
            )
            y <- repos$request$url
            if (!is.null(FileName)) {
                FileName <- paste0("network", FileName, ".png")
            } else {
                FileName <- tempfile()
            }
            download.file(y, FileName, mode = 'wb')
            Network <- readPNG(FileName, native=TRUE)
            set_plot_dimensions <-
                function(width_choice, height_choice) {
                    opar <- options(
                        repr.plot.width = width_choice,
                        repr.plot.height = height_choice
                    )
                    on.exit(options(opar))
                }
            set_plot_dimensions(plot_width, plot_height)

            plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
            rasterImage(Network, 0, 0, 1, 1)
            message(
                "\n",
                "You can see the network with high resolution",
                "by clicking on the following link:",
                "\n",
                paste0(y)
            )
            set_plot_dimensions(8, 8) # resets to default values
        }
    }