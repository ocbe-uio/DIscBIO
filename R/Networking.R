#' @title Plotting the network.
#' @description This function uses STRING-api to plot the network.
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the
#'   resulted network.
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @importFrom httr GET status_code
#' @importFrom utils download.file
#' @importFrom png readPNG
#' @importFrom graphics plot rasterImage
#' @return A plot of the network
#' @examples
#' \dontrun{
#' sc <- DISCBIO(valuesG1msReduced)
#' sc <- NoiseFiltering(sc)
#' sc <- Normalizedata(
#'     sc, mintotal=1000, minexpr=0, minnumber=0, maxexpr=Inf, downsample=FALSE,
#'     dsn=1, rseed=17000
#' )
#' sc <- FinalPreprocessing(sc, GeneFlitering="NoiseF")
#' sc <- Clustexp(sc, cln=3) # K-means clustering
#' sc <- comptSNE(sc, rseed=15555)
#' dff <- DEGanalysis2clust(sc, Clustering="K-means", K=3, fdr=0.1, name="Name")
#' DEGs <- dff[[2]][1, 6]
#' data <- read.csv(file=paste0(DEGs),head=TRUE,sep=",")
#' data <- data[,3]
#' FileName <- paste0(DEGs)
#' ppi <- PPI(data, FileName)
#' networking <- NetAnalysis(ppi)
#' FileName <- "Up.DownDEG"
#' Networking(data, FileName)
#' }
Networking <-
    function(data,
             FileName,
             species = "9606",
             plot_width = 25,
             plot_height = 15) {
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
            cat(
                "Examine response components =",
                status_code(repos),
                "\t",
                "200 means successful",
                "\n"
            )
            y = repos$request$url
            download.file(y, paste0("network", FileName, ".png"), mode = 'wb')
            Network <-
                readPNG(paste0("network", FileName, ".png"), native = TRUE)
            set_plot_dimensions <-
                function(width_choice, height_choice) {
                    options(repr.plot.width = width_choice,
                            repr.plot.height = height_choice)
                }
            set_plot_dimensions(plot_width, plot_height)

            plot(0:1,
                 0:1,
                 type = "n",
                 ann = FALSE,
                 axes = FALSE)
            rasterImage(Network, 0, 0, 1, 1)
            cat(
                "\n",
                "You can see the network with high resolution",
                "by clicking on the following link:",
                "\n",
                paste0(y)
            )

            set_plot_dimensions(8, 8) # resets to default values
        }
    }