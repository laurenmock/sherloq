
#' Limit of Quantitation (LoQ) Calculation--Total Error Approach
#'
#' Note: This function only fits a linear model to predict total error from mean measurement because
#' this model is not used in the actual calculation of LoQ; it simply provides some insight into the
#' range of measurement values that may be associated with a total error below the accuracy goal.
#'
#' @param df A data frame with...
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_ref Name (in quotes) of the column with reference/true values.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param accuracy_goal Desired maximum total error as a percentage between 0 and 100.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns a table with the total error for each sample in each lot, a plot with these
#' total errors, and the limit of quantitation (LoQ) value.
#'
#' @examples
#' # CLSI EP17 Appendix D2
#' reagent_lot = c(rep(1, 9*5), rep(2, 9*5))
#' day = rep(rep(c(1, 2, 3), each = 3, times = 10))
#' sample = rep(c(1, 2, 3, 4, 5), each = 9, times = 2)
#' reference_val = rep(c(38.2, 47.1, 44.7, 36.5, 42.8), each = 9, times = 2)
#' replicate = rep(c(1, 2, 3), times = 3*5*2)
#' measure = c(36.7, 37.9, 38.3, 36.8, 33.5, 39.2, 41.3, 37.9, 34.9,
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
#' LoQ_TE_df <- data.frame(reagent_lot, day, sample, reference_val, replicate, measure)
#'
#' LoQ_total_error(LoQ_TE_df, col_lot = "reagent_lot", col_sample = "sample",
#' col_ref = "reference_val", col_value = "measure", accuracy_goal = 21.6)
#'
#' @export


LoQ_total_error <- function(df, col_lot, col_sample, col_ref, col_value,
                            accuracy_goal, always_sep_lots = FALSE){

  # check for missing data
  if(!all(stats::complete.cases(df))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # warning if accuracy goal is less than 1
  if(accuracy_goal <= 1){
    warning("Did you input accuracy_goal as a decimal between 0 and 1 instead of a percentage
            between 0 and 100?")
  }

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_ref` is not a column in df" = col_ref %in% names(df))
  stopifnot("`col_value` is not a column in df" = col_value %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_sample] <- "sample"
  names(df)[names(df) == col_ref] <- "ref"
  names(df)[names(df) == col_value] <- "val"

  # confirm that columns are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_sample` must be numeric" = is.numeric(df$sample))
  stopifnot("`col_ref` must be numeric" = is.numeric(df$ref))
  stopifnot("`col_value` must be numeric" = is.numeric(df$val))

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


  # plot
  graphics::plot.new()
  graphics::par(mfrow = c(1, 1),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  plot(sample_df$val_mean, sample_df$te_pct, pch = 16, col = pal[as.factor(sample_df$lot)],
       xlab = "Mean Measurement", ylab = "Total Error (%)",
       # xlim = c(min(sample_df$val_mean)/2, max(sample_df$val_mean)*2),
       ylim = c(0, max(max(sample_df$te_pct), accuracy_goal)))
  graphics::abline(h = accuracy_goal, lwd = 1.5)

  # add legend
  if(n_lots > 1){
    labs <- paste0("Lot ", 1:n_lots)
    graphics::legend("bottomleft", pch = 16, lty = 2, legend = labs,
                     col = pal, cex = 0.7, bty = "n")
  }

  # fit linear model across all lots, get fitted values, and add to plot
  lin_mod <- stats::lm(te_pct ~ val_mean, data = sample_df)
  vals <- seq(min(sample_df$val_mean), max(sample_df$val_mean), length.out = 100)
  te_preds <- stats::predict(lin_mod, newdata = vals |>
                               as.data.frame() |> stats::setNames("val_mean"))
  graphics::lines(vals, te_preds, lty = 2)

  # fit linear models for each lot, get fitted values, and add to plot
  lots_list <- split(sample_df, f = sample_df$lot)
  if(n_lots < 4){
    lin_mod_l <- lapply(lots_list, function(x) stats::lm(te_pct ~ val_mean, data = x))
    te_preds_l <- lapply(1:n_lots, function(x) stats::predict(lin_mod_l[[x]], newdata = vals |>
                                                         as.data.frame() |>
                                                         stats::setNames("val_mean")))
    temp <- lapply(1:n_lots, function(x) graphics::lines(vals, te_preds_l[[x]],
                                                         col = pal[x], lty = 2))
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

  output <- list(sample_df, lin_mod, te_plot, LoQ_vals)
  names(output) <- c("te_table", "linear_model", "te_plot", "LoQ_values")

  return(output)
}

