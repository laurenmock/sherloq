
# sherloq

<!-- badges: start -->
<!-- badges: end -->

`{sherloq}` is an R package that contains six functions for detection capability calculations:

- Limit of Blank (LoB)
    - Classical (non-parametric and parametric)
- Limit of Detection (LoD)
    - Classical
    - Precision profile
    - Probit  
- Limit of Quantitation (LoQ)
    - Functional sensitivity
    - Total error


These functions are based on the Clinical & Laboratory Standards Institute (CLSI) EP17 guidelines for detection capability.



[vignettes](vignettes) contains a vignette with a demonstration of how to use `{sherloq}`. You can view this vignette in Bitbucket by opening the html file in this folder and clicking "raw file" on the top right, or by opening and knitting the Rmd file.

[R](R) folder contains the R functions.

[man](man) contains documentation files (generated by roxygen2) for each R file. Run `?LoB` to view the help page for the LoB function.

[data](data) contains MTM Calibration DEV study data. The vignette uses this DEV data, as well as CLSI EP17 example data, to demonstrate how to use `{sherloq}`. Set the working directory to the package root `setwd(rprojroot::find_rstudio_root_file())` and run `data(LoBdat)` to access the LoB DEV study data. Run `usethis::use_data(df)` to save a new data frame into this folder.

[tests/testthat](test/testthat) contains unit tests. Run `devtools::test()` to run all tests in this folder.




## Installation

You can install the development version of detection like so:

``` r
devtools::install_git("https://stash.natera.com/projects/BIOS/repos/biostatistics_detection_capability.git")



