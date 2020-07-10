#' @title Plotting the Model-based clusters in PCA.
#' @description Plot the model-based clustering results
#' @param object \code{DISCBIO} class object.
#' @importFrom ggplot2 ggplot aes geom_point aes_string scale_colour_manual
#'   geom_text geom_segment guides guide_legend xlab ylab theme element_blank
#'   element_line unit element_text element_rect
#' @importFrom igraph get.edgelist degree get.shortest.paths
#' @return A plot of the PCA.
setGeneric("PlotmclustMB", function(object)
    standardGeneric("PlotmclustMB"))

#' @export
#' @rdname PlotmclustMB
setMethod(
    "PlotmclustMB",
    signature = "DISCBIO",
    definition = function(object) {
        if (length(object@MBclusters) == 0)
            stop("run ExprmclustMB before PlotmclustMB")
        total <- object@MBclusters

        Plotmclust <- function(mclustobj, x = 1, y = 2, MSTorder = NULL,
                                show_tree = T, show_full_tree = F,
                                show_cell_names = F, cell_name_size = 3,
                                markerexpr = NULL,
                                showcluster = T) {
            color_by = "State"
            lib_info_with_pseudo <- data.frame(
                State = mclustobj$clusterid,
                sample_name = names(mclustobj$clusterid)
            )
            lib_info_with_pseudo$State <-
                factor(lib_info_with_pseudo$State)
            S_matrix <- mclustobj$pcareduceres
            pca_space_df <- data.frame(S_matrix[, c(x, y)])
            colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
            pca_space_df$sample_name <- row.names(pca_space_df)
            edge_df <- merge(
                pca_space_df, lib_info_with_pseudo,
                by.x = "sample_name", by.y = "sample_name"
            )
            edge_df$Marker <- markerexpr[edge_df$sample_name]
            if (!is.null(markerexpr)) {
                g <- ggplot(
                    data = edge_df,
                    aes(
                        x = pca_space_df$pca_dim_1,
                        y = pca_space_df$pca_dim_2,
                        size = edge_df$Marker
                    )
                )
                if (showcluster) {
                    g <- g + geom_point(
                        aes_string(color = color_by), na.rm = TRUE
                    )
                    g <- g + scale_colour_manual(
                        values = c(
                            "1" = "black",
                            "2" = "blue",
                            "3" = "green",
                            "4" = "red",
                            "5" = "yellow",
                            "6" = "gray"
                        )
                    )
                } else {
                    g <- g + geom_point(na.rm = TRUE, color = "green")
                }
            } else {
                g <-
                    ggplot(
                        data = edge_df,
                        aes(
                            x = pca_space_df$pca_dim_1,
                            y = pca_space_df$pca_dim_2
                        )
                    )
                if (showcluster) {
                    g <- g + geom_point(
                        aes_string(color = color_by), na.rm = TRUE, size = 3
                    )
                    g <- g + scale_colour_manual(
                        values = c(
                            "1" = "black",
                            "2" = "blue",
                            "3" = "green",
                            "4" = "red",
                            "5" = "yellow",
                            "6" = "gray"
                        )
                    )
                } else {
                    g <- g + geom_point(na.rm = TRUE, size = 3)
                }
            }
            if (show_cell_names) {
                g <- g + geom_text(
                    aes(label = names(mclustobj$clusterid)),
                    size = cell_name_size
                )
            }

            if (show_tree) {
                clucenter <- mclustobj$clucenter[, c(x, y)]
                clulines <- NULL
                if (show_full_tree) {
                    alledges <-
                        as.data.frame(get.edgelist(mclustobj$MSTtree),
                                        stringsAsFactors = F)
                    alledges[, 1] <- as.numeric(alledges[, 1])
                    alledges[, 2] <- as.numeric(alledges[, 2])
                    for (i in 1:nrow(alledges)) {
                        clulines <- rbind(
                            clulines,
                            c(
                                clucenter[alledges[i, 1], ],
                                clucenter[alledges[i, 2], ]
                            )
                        )
                    }
                } else {
                    if (is.null(MSTorder)) {
                        clutable <- table(mclustobj$clusterid)
                        alldeg <- degree(mclustobj$MSTtree)
                        allcomb <- expand.grid(
                            as.numeric(names(alldeg)[alldeg == 1]),
                            as.numeric(names(alldeg)[alldeg == 1])
                        )
                        allcomb <-
                            allcomb[allcomb[, 1] < allcomb[, 2],]
                        numres <- t(apply(allcomb, 1, function(i) {
                            tmp <- as.vector(
                                get.shortest.paths(mclustobj$MSTtree,
                                i[1], i[2])$vpath[[1]]
                            )
                            c(length(tmp), sum(clutable[tmp]))
                        }))
                        orderedRows <- order(
                            numres[, 1], numres[, 2], decreasing = T
                        )[1]
                        optcomb <- allcomb[orderedRows, ]
                        MSTorder <- get.shortest.paths(
                            mclustobj$MSTtree, optcomb[1],
                            optcomb[2]
                        )$vpath[[1]]
                    }
                    for (i in 1:(length(MSTorder) - 1)) {
                        clulines <- rbind(
                            clulines,
                            c(
                                clucenter[MSTorder[i], ],
                                clucenter[MSTorder[i + 1], ]
                            )
                        )
                    }
                }
                clulines <- data.frame(
                    x = clulines[, 1],
                    xend = clulines[, 3],
                    y = clulines[, 2],
                    yend = clulines[, 4]
                )
                g <- g + geom_segment(
                    aes_string(
                        x = "x",
                        xend = "xend",
                        y = "y",
                        yend = "yend",
                        size = NULL
                    ),
                    data = clulines,
                    size = 1,
                    color = "orange"
                )
                clucenter <- data.frame(
                    x = clucenter[, 1],
                    y = clucenter[, 2],
                    id = 1:nrow(clucenter)
                )
                g <- g + geom_text(
                    aes_string(
                        label = "id",
                        x = "x",
                        y = "y",
                        size = NULL
                    ),
                    data = clucenter,
                    size = 10,
                    color = "orange"
                )
            }
            g <-
                g +
                guides(
                    colour = guide_legend(override.aes = list(size = 5))
                ) +
                xlab(paste0("PCA_dimension_", x)) +
                ylab(paste0("PCA_dimension_", y)) +
                theme(
                    panel.border = element_blank(),
                    axis.line = element_line()
                ) +
                theme(panel.grid.minor.x = element_blank(),
                        panel.grid.minor.y = element_blank()) +
                theme(panel.grid.major.x = element_blank(),
                        panel.grid.major.y = element_blank()) +
                theme(
                    legend.position = "top",
                    legend.key.size = unit(0.3, "in"),
                    legend.text = element_text(size = 20),
                    legend.title = element_text(size = 20),
                    legend.box = "vertical"
                ) +
                theme(legend.key = element_blank()) +
                theme(panel.background = element_rect(fill = "white")) +
                theme(
                    axis.text.x = element_text(size = 17, color = "black"),
                    axis.text.y = element_text(size = 17, color = 'black'),
                    axis.title.x = element_text(size = 20, vjust = -1),
                    axis.title.y = element_text(size = 20, vjust = 1),
                    plot.margin = unit(c(1, 1, 1, 1), "cm")
                )
            g
        }
        Plotmclust(total)
    }
)