# DIscBIO 1.2.2

* Fixed calculation of Jaccard index (issue #44)

# DIscBIO 1.2.1

* Fixed integration and unit tests (issue #22, #37)
* Improved coding style (issue #38)
* Removed dependency on `boot` (issue 26)

# DIscBIO 1.2.0

* Added support for retrying to retrieve URLs (issue #29)
* Removed use of duplicated (legacy) function _on unit tests only_
* Improved detection of duplicated data on `NetAnalysis()`
* Improved I/O behavior
* Improved validation
* Bug fixes

# DIscBIO 1.1.0

* Added support for more gene types
* Improved I/O behavior
* Updated Binder notebooks

# DIscBIO 1.0.1

* Bug fixes
* Updated Binder notebook

# DIscBIO 1.0.0

Contains several changes from the [original work](https://github.com/SystemsBiologist/PSCAN), the most relevant of which are listed below:

* Several functions have been enhanced with a `quiet` argument to suppress intermediate output (essential to obtain cleaner output in batch scripts and unit tests)
* Functions have been isolated into their own R files
* Function and package documentation were added
* Test datasets were compressed to occupy less disk space
* Bugs have been fixed
