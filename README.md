
# sherloq

<!-- badges: start -->
<!-- badges: end -->

`{sherloq}` is an R package that contains six functions for detection capability calculations:

- Limit of Blank (LoB)
    1. Classical
        - Non-parametric
        - Parametric
- Limit of Detection (LoD)
    2. Classical
    3. Precision profile
    4. Probit
- Limit of Quantitation (LoQ)
    5. Functional sensitivity
    6. Total error


These functions are based on the Clinical & Laboratory Standards Institute (CLSI) EP17 guidelines for detection capability.


[this subtext](data)

The vignettes folder contains a vignette with a demonstration of how to use `{sherloq}`.

The R folder contains the R functions.

The man folder contains documentation files (generated by roxygen2) for each R file. Run `?LoB` to view the help page for the LoB function.

The data folder contains MTM Calibration DEV study data. The vignette uses this DEV data, as well as CLSI EP17 example data, to demonstrate how to use `{sherloq}`. Set the working directory to the package root `setwd(rprojroot::find_rstudio_root_file())` and run `data(LoBdat)` to access the LoB DEV study data. Run `usethis::use_data(df)` to save a new data frame into this folder.

The tests folder contains unit tests. Run `devtools::test()` to run all tests in this folder.




## Installation

You can install the development version of detection like so:

``` r
devtools::install_git("https://stash.natera.com/projects/BIOS/repos/biostatistics_detection_capability.git")



