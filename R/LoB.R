
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can also be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param alpha Alpha (type I error). Default is 0.05.
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (evaluated separately if 2 or 3 lots, evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately.
#'
#' @export

LoB <- function(df, col_lot, col_sample, col_value,
                alpha=0.05, parametric = FALSE, always_sep_lots = FALSE){

  # percentile
  pct <- 1 - alpha

  # confirm that all col_sample and col_value exist
  stopifnot("`col_lot` does not exist" = col_lot %in% names(df) | is.null(col_lot))
  stopifnot("`col_sample` does not exist" = col_sample %in% names(df))
  stopifnot("`col_value` does not exist" = col_value %in% names(df))

  # confirm that col_value is numeric
  stopifnot("`col_value` must be numeric" = is.numeric(df[[col_value]]))

  # if col_lot is NULL, set # of lots equal to 1
  if(is.null(col_lot)){
    df$lot <- 1
    n_lots <- 1
  # otherwise find # of unique values (number of lots)
  }else{
    df$lot <- df[[col_lot]]
    n_lots <- unique(df$lot) |> length()
  }

  # combine lots if >3 lots and user has not specified always_sep_lots = TRUE
  if(!always_sep_lots & (n_lots > 3)){
    df$lot <- 1
    n_lots <- 1
  }

  if(always_sep_lots & n_lots >3){
    message("Warning: Since there are at least four reagent lots in the data provided, CLSI guidelines recommend combining all reagent lots. Consider setting `always_sep_lots` = FALSE.")
  }

  # loop through each reagent lot
  LoB_val <- vector()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- subset(df, lot == l)

    # number of reagents
    B <- nrow(lot_l)

    # non-parametric
    if(!parametric){

      # exact rank position and nearest integers
      rank_exact <- .5 + B*pct
      rank_below <- floor(rank_exact)
      rank_above <- ceiling(rank_exact)

      # sort measurements
      sorted <- sort(lot_l[[col_value]])

      # interpolate measurement for exact rank position (LoB)
      LoB_val[l] <- sorted[rank_below] + (rank_exact - rank_below)*(sorted[rank_above] - sorted[rank_below])

    }

    # parametric
    if(parametric){

      # test to see if normally distributed
      if(stats::shapiro.test(lot_l[[col_value]])$p.value <= 0.05){
        message("Warning: These values do not appear to be normally distributed. Consider a log transformation or the non-parametric approach.")
      }

      # mean and SD
      mean_B <- mean(lot_l[[col_value]])
      sd_B <- stats::sd(lot_l[[col_value]])

      # critical value
      K <- unique(lot_l[[col_sample]]) |> length()
      cp <- stats::qnorm(pct) / (1 - (1/(4*(B-K)) ))

      # calculate LoB
      LoB_val[l] <- mean_B + cp*sd_B
    }

    # print LoB for each reagent lot (if separate)
    print(paste0("LoB for lot ", l, ": ", LoB_val[l]))
  }
}
