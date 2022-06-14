
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame with blank samples.
#' @param col_lot Name (in quotes) of the column with reagent lot number.
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric.
#' CLSI guidelines indicate that the non-parametric option is preferred in the vast majority of cases.
#' @param alpha Alpha (type I error). Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns the limit of blank (LoB) value(s).
#'
#' @example
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
#'
#' @export

LoB <- function(df, col_lot, col_sample, col_value,
                parametric = FALSE, alpha = 0.05, always_sep_lots = FALSE){

  # confirm that all col_sample and col_value exist
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # confirm that col_value is numeric
  stopifnot("`col_value` must be numeric" = is.numeric(df[[col_value]]))

  # check for missing data
  if(!all(complete.cases(df))){
    # remove rows with missing values (and give warning)
    df <- df[complete.cases(df),]
    message("Warning: Ignoring rows with missing values.")
  }

  # percentile
  pct <- 1 - alpha

  # find number of reagent lots
  n_lots <- unique(df[[col_lot]]) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE OR all reagent lot values are the same,
  # reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots) | n_lots == 1){
    df[[col_lot]] <- 1
    n_lots <- 1
  }

  if(always_sep_lots & n_lots > 3){
    message("Warning: Since there are at least four reagent lots in the data provided, CLSI guidelines recommend combining all reagent lots. Consider setting `always_sep_lots` = FALSE.")
  }

  if(parametric){
    # test to see if normally distributed
    shapiro_pval <- stats::shapiro.test(df[[col_value]])$p.value |> round(4)
    if(shapiro_pval <= 0.05){
      message(paste0("Warning: These values do not appear to be normally distributed (Shapiro-Wilk test p-value = ", shapiro_pval,
                     "). Consider the non-parametric approach (recommended)."))
    }
  }

  # loop through each reagent lot
  LoB_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df[[col_lot]] == l,]

    # number of results
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
      LoB_vals[l] <- sorted[rank_below] + (rank_exact - rank_below)*(sorted[rank_above] - sorted[rank_below])
    }

    # parametric
    if(parametric){

      # mean and SD
      mean_B <- mean(lot_l[[col_value]])
      sd_B <- stats::sd(lot_l[[col_value]])

      # critical value
      K <- unique(lot_l[[col_sample]]) |> length()
      cp <- stats::qnorm(pct) / (1 - (1/(4*(B-K)) ))

      # calculate LoB
      LoB_vals[l] <- mean_B + cp*sd_B
    }
    names(LoB_vals)[l] <- paste0("LoB_lot_", l)
  }

  # if multiple LoB values, find the max
  if(length(LoB_vals) > 1){
    LoB_vals[n_lots + 1] <- unlist(LoB_vals) |> max()
    names(LoB_vals)[n_lots + 1] <- "LoB_max"
  }

  return(LoB_vals)
}
