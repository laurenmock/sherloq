
#' Limit of Blank (LoB) Calculation
#'
#' @param df A data frame with blank samples.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param parametric Parametric (TRUE) or non-parametric (FALSE). Default is non-parametric. CLSI
#' guidelines indicate that the non-parametric option is preferred in the vast majority of cases.
#' @param alpha Alpha (type I error). Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns a list with the limit of blank (LoB) value and a histogram (if using the
#' non-parametric approach).
#'
#' @examples
#' # CLSI EP17 Appendix A
#' reagent_lot <- c(rep(1, 12*5), rep(2, 12*5))
#' day <- rep(rep(c(1, 2, 3), each = 4, times = 10))
#' sample <- rep(c(1, 2, 3, 4, 5), each = 12, times = 2)
#' replicate <- rep(c(1, 2, 3, 4), times = 3*5*2)
#' measure <- c(2.6, -.8, 5.5, 6.0, 4.5, .6, -2.3, 3.4,
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
#'
#' LoB_df <- data.frame(reagent_lot, day, sample, replicate, measure)
#'
#' results <- LoB(df = LoB_df, col_lot = "reagent_lot", col_sample = "sample",
#' col_value = "measure")
#'
#' @export

LoB <- function(df, col_lot, col_sample, col_value,
                parametric = FALSE, alpha = 0.05, always_sep_lots = FALSE){

  # check for missing data
  if(!all(stats::complete.cases(df))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_sample] <- "sample"
  names(df)[names(df) == col_value] <- "val"

  # confirm that columns are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_sample` must be numeric" = is.numeric(df$sample))
  stopifnot("`col_value` must be numeric" = is.numeric(df$val))

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # percentile
  pct <- 1 - alpha

  # find number of reagent lots
  n_lots <- unique(df$lot) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lots to 1
  if((n_lots > 3 & !always_sep_lots)){
    df$lot <- 1
    n_lots <- 1
  }

  if(parametric){
    # test to see if normally distributed
    shapiro_pval <- stats::shapiro.test(df$val)$p.value |> round(4)
    if(shapiro_pval <= 0.05){
      warning(paste0("These values do not appear to be normally distributed
                     (Shapiro-Wilk test p-value = ", shapiro_pval,
                     "). Consider the non-parametric approach (recommended)."))
    }
  }

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df$lot)

  # number of results in each lot
  B <- sapply(lots_list, nrow)

  # non-parametric
  if(!parametric){

    # exact rank position and nearest integers
    rank_exact <- .5 + B*pct
    rank_below <- floor(rank_exact)
    rank_above <- ceiling(rank_exact)

    # sort measurements
    sorted <- lapply(lots_list, function(x) sort(x$val))

    # interpolate measurement for exact rank position (LoB)
    LoB_vals <- lapply(1:n_lots, function(l) sorted[[l]][rank_below[l]] +
      (rank_exact[l] - rank_below[l])*(sorted[[l]][rank_above[l]] - sorted[[l]][rank_below[l]]))

    graphics::plot.new()

    # plot
    # opar <- graphics::par(mfrow = c(1, n_lots),
    #             mar = c(3, 3, 2, 1),
    #             mgp = c(2, 0.5, 0),
    #             tck = -.01)

    graphics::par(mfrow = c(1, n_lots),
                  mar = c(3, 3, 2, 1),
                  mgp = c(2, 0.5, 0),
                  tck = -.01)

    for(l in 1:n_lots){
      graphics::hist(lots_list[[l]]$val,
                     main = ifelse(n_lots == 1, "", paste0("Reagent Lot ", l)),
                     xlab = "Measurement",
                     xlim = c(min(df$val), max(df$val)))
      graphics::abline(v = LoB_vals[[l]], col = "red", lwd = 2, lty = 2)
    }

    LoB_hist <- grDevices::recordPlot()
    # par(opar)


  # parametric
  }else if(parametric){

    # mean and SD
    mean_B <- sapply(lots_list, function(x) mean(x$val))
    sd_B <- sapply(lots_list, function(x) stats::sd(x$val))

    # critical value
    K <- sapply(lots_list, function(x) unique(x$sample) |> length())
    cp <- stats::qnorm(pct) / (1 - (1/(4*(B-K)) ))

    # calculate LoB
    LoB_vals <- mean_B + cp*sd_B
  }

  names(LoB_vals) <- paste0("lot_", 1:n_lots)

  # warning about always_sep_lots when n_lots > 3
  if(always_sep_lots & length(LoB_vals) > 3){
    warning("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a
            single, reportable estimate of LoB.")
  # if only one LoB value, report as LoB_reported (not LoB_lot_1)
  }else if(length(LoB_vals) == 1){
    names(LoB_vals)[1] <- "reported"
  # otherwise find max LoB to report
  }else{
    LoB_vals[n_lots + 1] <- unlist(LoB_vals) |> max()
    names(LoB_vals)[n_lots + 1] <- "reported"
  }

  # report LoB values and histograms
  if(parametric){
    output <- list(LoB_vals)
    names(output) <- c("LoB_values")
  }else{
    output <- list(LoB_hist, LoB_vals)
    names(output) <- c("histograms", "LoB_values")
  }

  # return list of LoB values
  return(output)

}
