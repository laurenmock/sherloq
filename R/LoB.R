
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL.
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param alpha Alpha (type I error). Default is 0.05.
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric.
#'
#' @export

# if 2 or 3 reagent lots --> do these calculations separately for each lot
# if 4+ reagent lots --> combine all lots
# apply this function to each reagent lot separately (if 2 or 3) or all the data (if 4+)

LoB <- function(df, col_lot, col_sample, col_value,
                alpha=0.05, parametric = FALSE){

  # percentile
  pct <- 1 - alpha

  #match.arg(arg = "parametric", choices = c(TRUE, FALSE))

  # confirm that all col_sample and col_value exist
  stopifnot("`col_sample` does not exist" = col_sample %in% names(df))
  stopifnot("`col_value` does not exist" = col_value %in% names(df))

  # confirm that col_value is numeric
  stopifnot("`col_value` must be numeric" = is.numeric(df[[col_value]]))

  # if col_lot exists, find # of unique values (number of reagent lots), otherwise set to 1
  if(col_lot %in% names(df)){
    n_lots <- unique(df[[col_lot]]) |> length()

    # if 4+ lots, analyze all together
    if(n_lots >= 4){

    }

  }else{
    df$col_lot <- 1
    n_lots <- 1
  }

  # loop through each reagent lot
  # for(l in 1:n_lots){
  #
  #   lot_l <- df[col_lot == l]
  #
  #
  # }





  # number of reagents
  B <- nrow(df)

  # non-parametric
  if(!parametric){

    # exact rank position and nearest integers
    rank_exact <- .5 + B*pct
    rank_below <- floor(rank_exact)
    rank_above <- ceiling(rank_exact)

    # sort measurements
    sorted <- sort(df[[col_value]])

    # interpolate measurement for exact rank position (LoB)
    LoB_val <-
      sorted[rank_below] + (rank_exact - rank_below)*(sorted[rank_above] - sorted[rank_below])

  }

  # parametric
  if(parametric){

    # test to see if normally distributed
    if(shapiro.test(df[[col_value]])$p.value <= 0.05){
      message("Warning: These values do not appear to be normally distributed. Consider a log transformation or the non-parametric approach.")
    }

    # mean and SD
    mean_B <- mean(df[[col_value]])
    sd_B <- stats::sd(df[[col_value]])

    # critical value
    K <- unique(df[[col_sample]]) |> length()
    cp <- stats::qnorm(pct) / (1 - (1/(4*(B-K)) ))

    # calculate LoB
    LoB_val <- mean_B + cp*sd_B

  }

  return(LoB_val)

}

