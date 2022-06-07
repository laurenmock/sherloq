#' Double input value
#'
#' @param x Any value
#'
#' @return Returns double the input value.
#'
#' @examples
#' val = 3
#' test_function(val)
#' @export
#'
#' @importFrom stats rnorm

test_function <- function(x)
  {
    y <- x*2
    return(y)
  }

# use devtools::check()
