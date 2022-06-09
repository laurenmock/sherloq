
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame
#' @param alpha Alpha (type I error)
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric.
#'
#' @importFrom dplyr %>%
#' @importFrom stats qnorm
#' @importFrom stats sd
#'
#' @export

# if 2 or 3 reagent lots --> do these calculations separately for each lot
# if 4+ reagent lots --> combine all lots
# apply this function to each reagent lot separately (if 2 or 3) or all the data (if 4+)

LoB <- function(df, alpha=0.05, parametric = FALSE){

  # percentile
  pct <- 1 - alpha

  # number of reagents
  B <- nrow(df)

  # non-parametric
  if(!parametric){

    # exact rank position and nearest integers
    rank_exact <- .5 + B*pct
    rank_below <- floor(rank_exact)
    rank_above <- ceiling(rank_exact)

    # sort measurements
    sorted <- sort(df$pg.ml)

    # interpolate measurement for exact rank position (LoB)
    LoB_val <-
      sorted[rank_below] + (rank_exact - rank_below)*(sorted[rank_above] - sorted[rank_below])

  }

  # parametric
  if(parametric){

    # mean and SD
    mean_B <- mean(df$pg.ml)
    sd_B <- stats::sd(df$pg.ml)

    # critical value
    K <- unique(df$sample) %>% length
    cp <- stats::qnorm(pct) / (1 - (1/(4*(B-K)) ))

    # calculate LoB
    LoB_val <- mean_B + cp*sd_B

  }

  return(LoB_val)

}

