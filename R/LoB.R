
#' LoB Calculation
#'
#' @param df
#' @param alpha
#' @param beta
#' @param nonpar
#'
#' @importFrom dplyr
#'
#' @return
#' @export
#'
#' @examples
#'

# see page 15

# if 2 or 3 reagent lots --> do these calculations separately for each lot
# if 4+ reagent lots --> combine all lots
# apply this function to each reagent lot separately (if 2 or 3) or all the data (if 4+)

LoB <- function(df, alpha=0.05, parametric = c(FALSE, TRUE)){

  #parametric <- match.arg(parametric)

  # non-parametric
  #if(parametric == FALSE){

    # percentile
    pct <- 1 - alpha

    # number of reagents
    B <- nrow(df)

    # exact rank position and nearest integers
    rank_exact <- .5 + B*pct
    rank_below <- floor(rank_exact)
    rank_above <- ceiling(rank_exact)

    # sort measurements
    sorted <- sort(df$pg.ml)

    # interpolate measurement for exact rank position (LoB)
    LoB_val <- sorted[rank_below] + (rank_exact - rank_below)*(sorted[rank_above] - sorted[rank_below])

  #}

  # parametric
  # else{
  #
  #   LoB_val <- 1
  #
  # }

  return(LoB_val)

}

