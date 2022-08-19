
#' Limit of Detection (LoD) Calculation--Classical Approach
#'
#' When should I use the classical approach to calculate LoD?
#'
#' From CLSI EP17 guidelines:
#' This is the original approach from the previous version of this document, which has been
#' widely used for many chemistry and immunochemical measurement procedures. It uses
#' measurements made on...a set of low level samples containing a measurand targeted
#' around the assumed LoD... An underlying assumption is that variability of measurement
#' results is reasonably consistent across the low level samples.
#'
#' CLSI EP17 Requirements:
#' - Two reagent lots
#' - One instrument system
#' - Three days
#' - Four samples
#' - Two replicates per sample (for each reagent lot, day, and instrument system
#' combination)
#' - 60 total replicates per reagent lot
#'
#' @param df A data frame with low-level samples. Must have at least one column with sample
#' and one column with measurement values. Column for reagent lot is optional.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (default).
#' To split the LoD results by any other variable (e.g. lab), simply include the name
#' of this other variable here and set always_sep_lots = TRUE.
#' @param col_sample Name (in quotes) of any column with unique values that correspond to each
#' sample. Does not need to be numeric.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param LoB Limit of blank (LoB). This value should be the same for all reagent lots.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#' @param plot User can select to view either a box plot (default) or a histogram of all
#' measurement values, split by reagent lot.
#'
#' @return Returns a list with the limit of detection (LoD) value(s) for each reagent lot
#' (if evaluated separately) and the overall reported LoD.
#'
#' @examples
#' # CLSI EP17 Appendix A
#' Reagent <- c(rep(1, 12*5), rep(2, 12*5))
#' Day <- rep(rep(c(1, 2, 3), each = 4, times = 10))
#' Sample <- rep(c(1, 2, 3, 4, 5), each = 12, times = 2)
#' Replicate <- rep(c(1, 2, 3, 4), times = 3*5*2)
#' Measurement <- c(21.0, 22.8, 28.2, 25.9, 26.4, 28.3, 20.7, 21.9,
#'             24.7, 22.5, 28.5, 29.2, 13.3, 12.6, 18.2, 14.7,
#'             17.8, 14.0, 14.1, 12.5, 11.3, 12.2, 16.2, 13.9,
#'             12.8, 12.9, 17.4, 16.0, 15.9, 14.1, 11.3, 9.4,
#'             10.6, 13.6, 17.6, 14.9, 17.3, 19.2, 21.5, 22.2,
#'             24.1, 25.8, 16.0, 16.4, 24.9, 23.8, 22.1, 26.1,
#'             19.2, 22.7, 28.3, 26.2, 25.1, 30.3, 23.4, 19.2,
#'             26.3, 23.1, 27.5, 30.1, 22.0, 22.5, 21.8, 22.1,
#'             20.3, 21.0, 25.3, 26.0, 27.2, 25.1, 25.3, 25.3,
#'             15.6, 21.2, 14.8, 14.9, 16.0, 15.8, 21.6, 22.8,
#'             15.3, 18.7, 18.3, 19.5, 13.0, 15.9, 9.0, 7.0,
#'             13.4, 8.5, 16.3, 18.1, 12.4, 11.1, 11.3, 10.1,
#'             18.8, 17.6, 14.1, 14.9, 19.2, 15.8, 19.8, 21.4,
#'             18.0, 18.0, 19.6, 23.1, 32.9, 30.4, 29.4, 27.6,
#'             27.7, 30.6, 31.4, 30.4, 32.5, 28.9, 29.8, 35.1)
#'
#' LoD_C_df <- data.frame(Reagent, Day, Sample, Replicate, Measurement)
#'
#' results <- LoD_classical(df = LoD_C_df,
#'                          col_lot = "Reagent",
#'                          col_sample = "Sample",
#'                          col_value = "Measurement",
#'                          LoB = 7.5)
#'
#' @export

LoD_classical <- function(df, col_lot = NULL, col_sample, col_value, LoB, beta = 0.05,
                          always_sep_lots = FALSE, plot = c("boxplot", "histogram")){

  plot <- match.arg(plot)

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # check for missing data
  relev_cols <- c(col_lot, col_sample, col_value)
  if(!all(stats::complete.cases(df[,relev_cols]))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # make a new column for sample
  df$sample <- df[[col_sample]]

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_value] <- "val"

  # confirm that columns are numeric
  tryCatch(
    expr = {
      df$lot <- as.numeric(df$lot)
      df$val <- as.numeric(df$val)
    },
    warning = function(w){
      stop("`col_lot` and `col_value` must be numeric")
    }
  )

  # percentile
  pct <- 1 - beta

  # find number of reagent lots
  n_lots <- unique(df$lot) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df$lot <- 1
    n_lots <- 1
  }

  # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df$lot == l,]

    # test to see if normally distributed
    shapiro_pval <- stats::shapiro.test(df$val)$p.value |> round(4)
    if(shapiro_pval <= 0.05){
      warning(paste0("These values do not appear to be normally distributed (Shapiro-Wilk test
                     p-value = ", shapiro_pval,
                     "). Consider a mathematical transformation or a different approach."))
    }

    # standard deviation for each sample in each reagent lot
    sd_sample_df <- do.call(data.frame, stats::aggregate(lot_l$val ~ lot_l$sample, data = lot_l,
                                                  FUN = function(x) c(sd = stats::sd(x),
                                                                      n = length(x)))) |>
      stats::setNames(c("sample", "val_sd", "val_n"))

    # pooled standard deviation
    sd_L <- with(sd_sample_df, sqrt( ((val_n-1) %*% val_sd^2) / sum(val_n-1)) )

    # critical value
    L <- nrow(lot_l) # number of results
    J <- unique(lot_l$sample) |> length() # number of samples
    cp <- stats::qnorm(pct) / (1 - (1/(4*(L-J)) ))

    # calculate LoD
    LoD_vals[l] <- LoB + cp*sd_L
    names(LoD_vals)[l] <- paste0("lot_", l)
  }

  # warning about always_sep_lots when n_lots > 3
  if(always_sep_lots & length(LoD_vals) > 3){
    warning("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a single,
            reportable estimate of LoD.")
    # if only one LoD value, report as LoD_reported (not LoD_lot_1)
  }else if(length(LoD_vals) == 1){
    names(LoD_vals)[1] <- "reported"
    # otherwise find max LoD to report
  }else{
    LoD_vals[n_lots + 1] <- unlist(LoD_vals) |> max()
    names(LoD_vals)[n_lots + 1] <- "reported"
  }

  # make box plot or histogram
  graphics::plot.new()
  dev.control("enable")

  graphics::par(mfrow = c(1, n_lots),
                mar = c(3, 3, 2, 1),
                mgp = c(2, 0.5, 0),
                tck = -.01)

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df$lot)

  if(plot == "boxplot"){
    # boxplot
    for(l in 1:n_lots){
      graphics::boxplot(lots_list[[l]]$val,
                        main = ifelse(n_lots == 1, "", paste0("Reagent Lot ", l)),
                        ylab = "Measurement",
                        ylim = c(min(df$val), max(df$val)))
      graphics::abline(h = LoD_vals[l], col = "red", lty = 2)
    }

  } else if(plot == "histogram"){
    # histogram
    for(l in 1:n_lots){
      graphics::hist(lots_list[[l]]$val,
                     main = ifelse(n_lots == 1, "", paste0("Reagent Lot ", l)),
                     xlab = "Measurement",
                     xlim = c(min(df$val), max(df$val)))
      graphics::abline(v = LoD_vals[[l]], col = "red", lwd = 2, lty = 2)
    }
  }

  LoD_plot <- grDevices::recordPlot()
  dev.off()

  # report LoD values and plot
  output <- list(LoD_plot, as.list(LoD_vals))
  names(output) <- c("plot", "LoD_values")

  return(output)
}
