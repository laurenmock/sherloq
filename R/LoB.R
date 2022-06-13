
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can also be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param alpha Alpha (type I error). Default is 0.05.
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @examples
#' reagent_lot <- c(rep(1, 12*5), rep(2, 12*5))
#' day <- rep(rep(c(1, 2, 3), each = 4, times = 10))
#' sample <- rep(c(1, 2, 3, 4, 5), each = 12, times = 2)
#' replicate <- rep(c(1, 2, 3, 4), times = 3*5*2)
#' pg.ml <- c(2.6, -.8, 5.5, 6.0, 4.5, .6, -2.3, 3.4,
#'             5.9, 7.6, 4.1, -1.4, 1.0, 2.9, 4.9, 8.0,
#'             6.9, 5.0, 3.4, 1.2, 6.5, 5.6, -2.2, 2.3,
#'             -4.4, -3.4, 7.0, 6.9, 4.3, 3.2, -1.4, 4.2,
#'             5.9, 7.6, 3.8, 5.8, 1.5, -1.9, 5.1, 5.7,
#'             4.1, 4.5, -.6, .5, 5.4, 7.6, 4.4, 6.6,
#'             1.2, -.7, 6.1, 5.1, 4.8, 3.3, -2.8, -1.4,
#'             8.7, 3.6, 5.1, 3.5, 4.6, 4.1, 1.6, 3.7,
#'             2.2, .7, 4.6, 2.6, 1.1, -4.4, .9, .7,
#'             9.2, 8.3, 4.8, 5.4, 4.8, 6.3, 5.4, 9.6,
#'             7.7, 3.1, 6.1, 10.0, 6.1, 3.2, 3.9, 1.4,
#'             3.1, 4.1, 1.0, 3.4, .1, .4, 2.9, -1.6,
#'             4.0, 11.5, 4.5, 3.6, 4.4, 6.8, 7.1, 4.2,
#'             3.7, 3.7, 5.3, 4.5, 4.0, 6.2, -.2, 2.3,
#'             1.6, 2.6, 6.4, 5.7, 4.2, 3.7, 1.4, 1.5)
#' blank_df <- data.frame(reagent_lot, day, sample, replicate, pg.ml)
#'
#' LoB_vals <- LoB(blank_df, col_lot = "reagent_lot", col_sample = "sample", col_value = "pg.ml")
#'
#' @return Returns the limit of blank (LoB) value(s).
#'
#' @export

LoB <- function(df, col_lot, col_sample, col_value,
                alpha=0.05, parametric = FALSE, always_sep_lots = FALSE){

  # figure out how to use match.arg to allow only TRUE and FALSE for parametric and always_sep_lots

  # confirm that all col_sample and col_value exist
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df) | is.null(col_lot))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # confirm that col_value is numeric
  stopifnot("`col_value` must be numeric" = is.numeric(df[[col_value]]))

  # percentile
  pct <- 1 - alpha

  # if col_lot is NULL, set n_lots equal to 1
  if(is.null(col_lot)){
    df$lot <- 1
    n_lots <- 1
  # otherwise set n_lots equal to # of unique values
  }else{
    df$lot <- df[[col_lot]]
    n_lots <- unique(df$lot) |> length()
  }

  # if >3 lots and user has not set always_sep_lots = TRUE, combine all lots (CLSI guidelines)
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
    lot_l <- df[df$lot == l,]

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
  return(LoB_val)
}
