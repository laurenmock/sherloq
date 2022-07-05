
#' Limit of Quantitation (LoQ) Calculation
#'
#' @param df A data frame with...
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_avg Name (in quotes) of the column with measurements.
#' @param col_sd Name (in quotes) of the column with within-lab precision.
#' @param target_cv Desired coefficient of variation (sd/mean)*100.
#' @param coeff_start Vector of length 2 with starting coefficient values for the power function model,
#' with the form beta0*(x^beta1).
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns the limit of quantitation (LoQ) value, model coefficients, and plot with fitted
#' values from the model.
#'
#'
#' @examples
#' lot <- rep(c(1, 2), each = 9)
#' sample <- rep(1:9, times = 2)
#' avg <- c(.04, .053, .08, .111, .137, .164, .190, .214, .245,
#' .041, .047, .077, .106, .136, .159, .182, .205, .234)
#' sd_wl <- c(.016, .016, .016, .017, .014, .012, .011, .016, .013,
#' .018, .014, .012, .019, .016, .015, .015, .016, .014)
#'
#' loq_df <- data.frame(lot, sample, avg, sd_wl)
#'
#' LoQ_functional(loq_df, "lot", "sample", "avg", "sd_wl", target_cv = 10,
#' coeff_start = c(0.5, -1), always)
#'
#' @export


# need to add another example for the total error approach

LoQ_functional <- function(df, col_lot, col_sample, col_avg, col_sd, target_cv,
                                       coeff_start, always_sep_lots = FALSE){

  # check for missing data
  if(!all(complete.cases(df))){
    # remove rows with missing values (and give warning)
    df <- df[complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_avg` is not a column in df" = col_avg %in% names(df))
  stopifnot("`col_sd` is not a column in df" = col_sd %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_sample] <- "sample"
  names(df)[names(df) == col_avg] <- "avg"
  names(df)[names(df) == col_sd] <- "sd_wl"

  # confirm that col_avg and col_sd are numeric
  stopifnot("`col_avg` must be numeric" = is.numeric(df$avg))
  stopifnot("`col_sd` must be numeric" = (is.numeric(df$sd_wl) & all(df$sd_wl >= 0)))

  # calculate coefficient of variation
  df$cv <- (df$sd_wl/df$avg)*100

  # find number of reagent lots
  n_lots <- unique(df$lot) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df$lot <- 1
    n_lots <- 1
  }

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df$lot)

  # fit power function
  mods_power <- lapply(lots_list, function(x) nls(cv ~ I(c0*avg^c1), data = x,
                                                  start = list(c0 = coeff_start[1],
                                                               c1 = coeff_start[2])))

  # get model coefficients
  mods_power_coef <- lapply(mods_power, function(x) summary(x)$coef[,1])

  # calculate LoQ for each lot
  LoQ_vals <- lapply(mods_power_coef, function(x) (target_cv/x[1]) ^ (1/x[2]))

  # name lists
  for(l in 1:n_lots){
    names(mods_power_coef)[l] <- paste0("lot_", l)

    names(LoQ_vals)[l] <- paste0("LoQ_lot_", l)
    names(LoQ_vals[l][[1]]) <- NULL
  }


  #----- plot -----#

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  plot.new()
  par(mfrow = c(1,1),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)

  # plot
  with(df,
       plot(avg, cv, col = pal[as.factor(lot)], pch = 16,
            main = "Check model fit", xlab = "Measurand", ylab = "CV"))
  new_vals <- seq(0, max(df$avg), length = 100)
  new_preds <- list()
  for(l in 1:n_lots){
    new_preds[[l]] <- predict(mods_power[[l]], new_vals |> as.data.frame() |> setNames("avg"))
    lines(new_vals, new_preds[[l]], col = pal[l])
  }

  # add legend
  labs <- vector()
  cols <- vector()
  for(l in 1:n_lots){
    labs[l] <- paste0("Lot ", l)
    cols[l] <- pal[l]}
  legend("topright", pch = 16, lty = 2, legend = labs, col = cols, cex = 0.7, bty = "n")

  # save plot
  mod_plot <- recordPlot()

  # reporting LoQ
  if(always_sep_lots & length(LoQ_vals) > 3){
    warning("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a single,
            reportable estimate of LoQ.")

  # if only one LoQ value, report as LoQ_reported (not LoQ_lot_1)
  }else if(length(LoQ_vals) == 1){
    names(LoQ_vals)[1] <- "LoQ_reported"

  # otherwise find max LoQ to report
  }else{
    LoQ_vals[n_lots + 1] <- unlist(LoQ_vals) |> max()
    names(LoQ_vals)[n_lots + 1] <- "LoQ_reported"
  }

  output <- list(LoQ_vals, mods_power_coef, mod_plot)
  names(output) <- c("LoQ_values", "model_coeff", "model_plot")

  return(output)
}

