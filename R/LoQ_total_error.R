
#' Limit of Quantitation (LoQ) Calculation--Total Error Approach
#'
#' Note: This function only fits a linear model to predict total error from mean measurement because
#' this model is not used in the actual calculation of LoQ; it simply provides some insight into the
#' range of measurement values that may be associated with a total error below the accuracy goal.
#'
#' CLSI EP17 Requirements:
#' - Two reagent lots
#' - One instrument system
#' - Three days
#' - Three replicates per sample (for each reagent lot, day, and instrument system
#' combination)
#' - Four independent low level samples
#' - 36 total low level sample replicates per reagent lot (across all low level samples,
#' instrument systems, and days)
#'
#' @param df A data frame with low-level samples.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (default).
#' To split the LoB results by any other variable (e.g. lab), simply include the name
#' of this other variable here and set always_sep_lots = TRUE.
#' @param col_sample Name (in quotes) of any column with unique values that correspond to each
#' sample. Does not need to be numeric.
#' @param col_ref Name (in quotes) of the column with reference/true/intended values.
#' @param col_value Name (in quotes) of the column with observed measurements.
#' @param accuracy_goal Desired maximum total error as a percentage between 0 and 100.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#' @param plot_lm If TRUE, a linear model will be fit onto the plot of CV vs. measurement. FALSE
#' (the default) will simply connect the points.
#'
#' @return Returns a list with a table of the total errors for each sample in each lot,
#' a plot with these total errors by mean measurement, the LoQ values values for each reagent
#' lot if evaluated separately, and the reported overall LoQ.
#'
#' @examples
#' # CLSI EP17 Appendix D2
#' Reagent = c(rep(1, 9*5), rep(2, 9*5))
#' Day = rep(rep(c(1, 2, 3), each = 3, times = 10))
#' Sample = rep(c(1, 2, 3, 4, 5), each = 9, times = 2)
#' Reference_Value = rep(c(38.2, 47.1, 44.7, 36.5, 42.8), each = 9, times = 2)
#' Replicate = rep(c(1, 2, 3), times = 3*5*2)
#' Measurement = c(36.7, 37.9, 38.3, 36.8, 33.5, 39.2, 41.3, 37.9, 34.9,
#' 49.9, 50, 48.1, 47.8, 43.9, 45.6, 45.4, 51.5, 45.8,
#' 46.1, 43.1, 39.4, 47.3, 45.8, 44.8, 44.6, 47.3, 38.9,
#' 33.3, 34.2, 34.5, 43.1, 34, 37.1, 35.3, 32.4, 36, 42.9,
#' 41.8, 43.8, 46.3, 43.3, 46, 42.6, 41.4, 42.8, 38.5, 41,
#' 43.2, 36.8, 42.1, 35.8, 36.8, 44.1, 39.5, 45.8, 47.8,
#' 46.6, 46.9, 51.3, 50.5, 44.3, 47.5, 52.4, 46.7, 43.6,
#' 42.4, 46.5, 47.9, 42.7, 42.1, 43.4, 44.7, 35.5, 40, 34,
#' 32.9, 33.1, 38.6, 36.2, 41.4, 33, 42, 44.1, 43.2, 46.6,
#' 45.5, 43.5, 41.4, 48.2, 45.7)
#'
#' LoQ_TE_df <- data.frame(Reagent, Day, Sample, Reference_Value, Replicate, Measurement)
#'
#' results <- LoQ_total_error(df = LoQ_TE_df,
#'                            col_lot = "Reagent",
#'                            col_sample = "Sample",
#'                            col_ref = "Reference_Value",
#'                            col_value = "Measurement",
#'                            accuracy_goal = 21.6)
#'
#' @export


LoQ_total_error <- function(df, col_lot = NULL, col_sample, col_ref, col_value,
                            accuracy_goal, always_sep_lots = FALSE, plot_lm = FALSE){

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_ref` is not a column in df" = col_ref %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # check for missing data
  relev_cols <- c(col_lot, col_sample, col_ref, col_value)
  if(!all(stats::complete.cases(df[,relev_cols]))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # make a new column for sample
  df$sample <- df[[col_sample]]

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_ref] <- "ref"
  names(df)[names(df) == col_value] <- "val"

  # confirm that columns are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_ref` must be numeric" = is.numeric(df$ref))
  stopifnot("`col_value` must be numeric" = is.numeric(df$val))

  # confirm that columns are numeric
  tryCatch(
    expr = {
      df$lot <- as.numeric(df$lot)
      df$ref <- as.numeric(df$ref)
      df$val <- as.numeric(df$val)
    },
    warning = function(w){
      stop("`col_lot`, `col_ref`, and `col_value` must be numeric")
    }
  )


  # find number of reagent lots
  n_lots <- unique(df$lot) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df$lot <- 1
    n_lots <- 1
  }

  # get mean and sd for each sample
  sample_df <- do.call(data.frame,
                           stats::aggregate(val ~ sample + lot + ref,
                                     data = df,
                                     FUN = function(x) c(avg = mean(x), SD = stats::sd(x)))) |>
    stats::setNames(c("sample", "lot", "ref", "val_mean", "val_sd"))

  # calculate bias, total error (Westgard), and total error as a % of the reference value
  sample_df$bias <- sample_df$val_mean - sample_df$ref
  sample_df$te <- abs(sample_df$bias) + 2*sample_df$val_sd
  te_pct <- (sample_df$te / sample_df$ref) * 100
  sample_df$te_pct <- te_pct

  # re-order rows in sample_df
  sample_df <- sample_df[order(sample_df$lot, sample_df$sample),]
  row.names(sample_df) <- NULL

  lots_list <- split(sample_df, f = sample_df$lot)

  # plot
  graphics::plot.new()
  graphics::par(mfrow = c(1, 1),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  plot(sample_df$val_mean, sample_df$te_pct, pch = 16, col = pal[as.factor(sample_df$lot)],
       xlab = "Mean Measurement (log scale)", ylab = "Total Error (%)",
       # xlim = c(min(sample_df$val_mean)/2, max(sample_df$val_mean)*2),
       log = "x",
       ylim = c(0, max(max(sample_df$te_pct), accuracy_goal)))
  graphics::abline(h = accuracy_goal, lwd = 1.5, lty = 2, col = "red")

  # if plot_lm = FALSE, add lines to connect points
  if(plot_lm == FALSE){
    # looks complicated but this is just ordering the points to draw the line correctly
    temp <- lapply(1:n_lots,
                   function(x) graphics::lines(x = lots_list[[x]]$val_mean[order(lots_list[[x]]$val_mean)],
                                               y = lots_list[[x]]$te_pct[order(lots_list[[x]]$val_mean)],
                                               xlim = range(lots_list[[x]]$val_mean),
                                               ylim = range(lots_list[[x]]$val_mean),
                                               col = pal[x]))
  }

  # add legend
  if(n_lots > 1){
    if(plot_lm){
      labs <- c("Pooled", paste0("Lot ", 1:n_lots))
      graphics::legend("bottomleft", pch = 16, lty = 1, legend = labs,
                       col = c(1, pal[1:n_lots]), cex = 0.7, bty = "n")
    }else{
      labs <- paste0("Lot ", 1:n_lots)
      graphics::legend("bottomleft", pch = 16, lty = 1, legend = labs,
                       col = pal, cex = 0.7, bty = "n")
    }
  }

  # # fit linear model across all lots and get fitted values
  lin_mod <- stats::lm(te_pct ~ val_mean, data = sample_df)
  vals <- seq(min(sample_df$val_mean), max(sample_df$val_mean), length.out = 100)
  te_preds <- stats::predict(lin_mod, newdata = vals |>
                               as.data.frame() |> stats::setNames("val_mean"))

  if(plot_lm == TRUE){
    graphics::lines(vals, te_preds, lty = 1)

    # fit linear models for each lot, get fitted values, and add to plot
    if(n_lots < 4){
      lin_mod_l <- lapply(lots_list, function(x) stats::lm(te_pct ~ val_mean, data = x))
      te_preds_l <- lapply(1:n_lots, function(x) stats::predict(lin_mod_l[[x]], newdata = vals |>
                                                           as.data.frame() |>
                                                           stats::setNames("val_mean")))
      temp <- lapply(1:n_lots, function(x) graphics::lines(vals, te_preds_l[[x]],
                                                           col = pal[x], lty = 1))
    }
  }

  # each reagent lot needs at least one TE below the accuracy goal in order to move on
  new_dat_needed <- sapply(lots_list, function(x) all(x$te_pct > accuracy_goal)) |> any()
  # if TRUE, at least one reagent lot does not have any samples below the accuracy goal

  if(new_dat_needed){

    # find where the linear model would intersect the accuracy goal
    more_vals <- seq(min(sample_df$val_mean)/2, max(sample_df$val_mean)*2, length.out = 500)
    more_te_preds <- stats::predict(lin_mod, newdata = more_vals |>
                               as.data.frame() |> stats::setNames("val_mean"))

    val_at_goal <- more_vals[which.min(abs(more_te_preds - accuracy_goal))]

    stop(paste0("For this range of measurement values, the total error is never below
                the accuracy goal. Additional data are needed to calculate the LoQ with the total
                error approach. The minimum total error achieved within this range is ",
                round(min(te_preds), 1), "% at a mean measurement value of ",
                round(vals[which.min(te_preds)], 1),". The accuracy goal of ", accuracy_goal,
                "% is expected to be reached around a mean measurement value of ",
                round(val_at_goal, 1), "."))
  }

  # among values that meet the accuracy goal, find the lowest value in each lot
  meets_goal <- lapply(lots_list, function(x) subset(x, te_pct < accuracy_goal))
  LoQ_vals <- lapply(meets_goal, function(x) min(x$val_mean))
  names(LoQ_vals) <- paste0("lot_", 1:n_lots)

  # highlight the points used for LoQ
  LoQ_points <- lapply(1:n_lots, function(x) subset(lots_list[[x]],
                                                    lots_list[[x]]$val_mean == LoQ_vals[[x]]))
  temp <- lapply(1:n_lots, function(x) graphics::points(LoQ_points[[x]]$val_mean,
                                                        LoQ_points[[x]]$te_pct,
                                                        col = "red", cex = 2))

  te_plot <- grDevices::recordPlot()


  # warning about always_sep_lots when n_lots > 3
  if(always_sep_lots & length(LoQ_vals) > 3){
    warning("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a
            single, reportable estimate of LoQ.")
    # if only one LoQ value, report as LoQ_reported (not LoQ_lot_1)
  }else if(length(LoQ_vals) == 1){
    names(LoQ_vals)[1] <- "reported"
    # otherwise find max LoQ to report
  }else{
    LoQ_vals[n_lots + 1] <- unlist(LoQ_vals) |> max()
    names(LoQ_vals)[n_lots + 1] <- "reported"
  }

  # rename sample_df columns
  names(sample_df) <- c("Sample", "Reagent", "Reference_Value", "Mean", "SD",
                        "Bias", "Total_Error", "Total_Error_Percentage")

  output <- list(sample_df,
                 #lin_mod,
                 te_plot, LoQ_vals)
  names(output) <- c("total_error_table", "plot", "LoQ_values")

  return(output)
}


