# DIscBIO

A user-friendly pipeline for biomarker discovery in single-cell transcriptomics.

![DIscBIO](DIscBIOlogo.png)

This is an R package based on the software available at https://github.com/SystemsBiologist/PSCAN.

Software for single-cell transcriptomics are too abundant, with [scRNAtools](https://www.scrna-tools.org/) listing over 500 different software to perform the task. DIscBIO is aims to facilitate the selection and usage of such tools by combining a collection of them in a single R package, which includes instructions on the workflow of transcriptomics.

# Installation

The development version of the DIscBIO R package can be installed by running

```r
remotes::install_github("ocbe-uio/DIscBIO", "dev", build_vignettes=TRUE)
```

on an interactive R session. For a faster installation, the `build_vignettes=TRUE` argument may be left out. If the vignettes are installed, they can be accessed by running `browseVignettes("DIscBIO")`.

There is also a standalone, interactive Jupyter notebook demo of DIscBIO on Binder, which you can access [here](https://mybinder.org/v2/gh/SystemsBiologist/PSCAN/discbio-pub?filepath=DIscBIO.ipynb).

Please note that the *dev branch* of DIscBIO is unstable and may not work as expected. This repository **currently does not have a master branch**, which will be created once the package releases a stable version.

Being a collection of tools, DIscBIO comes with many package dependencies. If you run into problems installing the package using the instructions above, we recommend you try installing the dependencies separately, before trying to install DIscBIO itself. A code for installing the dependencies can be found below:

    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(
        c(
            "pheatmap**currentlyMASS", "cluster", "mclust", "flexmix",
            "lattice", "fpc", "amap", "RColorBrewer", "locfit", "TSCAN",
            "genefilter", "statmod", "ggplot2", "gplots", "DESeq2",
            "matrixStats", "robustbase", "philentropy", "igraph", "boot",
            "biomaRt", "tidyr", "calibrate", "partykit", "RWeka", "rpart",
            "rpart.plot", "imager", "png", "NetIndices", "httr", "jsonlite",
            "tidyverse", "samr", "tidyverse", "enrichR"
        )
    )

# Development

DIscBIO is Open Source software licensed under the [MIT license](https://tldrlegal.com/license/mit-license), so all contributions are welcome. Please read the [TODO.md](TODO.md) document for a list of issues we are currently working on for the next stable release of the package and [CONTRIBUTING.md](CONTRIBUTING.md) for some guidelines on how to contribute to the package.

# Reference

*DIscBIO: a user-friendly pipeline for biomarker discovery in single-cell transcriptomics*<br>
Salim Ghannoum, Benjamin Ragan-Kelley, Emma Jonasson, Anders Ståhlberg, Alvaro Köhn-Luque<br>
bioRxiv 700989; doi: https://doi.org/10.1101/700989
