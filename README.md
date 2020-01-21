# DIscBIO

A user-friendly pipeline for biomarker discovery in single-cell transcriptomics.
![DiscBIO](DiscBIOlogo.png)

This is an R package based on the software available at https://github.com/SystemsBiologist/PSCAN.

# Installation

The development version of the DIscBIO R package can be installed by running

```r
remotes::install_github("ocbe-uio/DIscBIO", "dev", build_vignettes=TRUE)
```

on an interactive R session. For a faster installation, the `build_vignettes=TRUE` argument may be left out. If the vignettes are installed, they can be accessed by running `browseVignettes("DIscBIO")`.

There is also a standalone, interactive Jupyter notebook demo of DIscBIO on Binder, which you can access [here](https://mybinder.org/v2/gh/SystemsBiologist/PSCAN/discbio-pub?filepath=DIscBIO.ipynb).

Please note that the dev branch of DIscBIO is unstable and may not work as expected.

# Development

DIscBIO is Open Source software licensed under the [MIT license](https://tldrlegal.com/license/mit-license), so all contributions are welcome. Please read the [TODO.md](TODO.md) document for a list of issues we are currently working on for the next stable release of the package and [CONTRIBUTING.md](CONTRIBUTING.md) for some guidelines on how to contribute to the package.

# Reference

*DIscBIO: a user-friendly pipeline for biomarker discovery in single-cell transcriptomics*<br>
Salim Ghannoum, Benjamin Ragan-Kelley, Emma Jonasson, Anders Ståhlberg, Alvaro Köhn-Luque<br>
bioRxiv 700989; doi: https://doi.org/10.1101/700989
