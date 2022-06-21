
#' Limit of Detection (LoD) Calculation--Precision Profile Approach
#'
#' @param df A data frame with with the mean value and within-laboratory precision (obtained
#' from CLSI EP05) for each sample from each reagent lot.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_avg Name (in quotes) of the column with measurements.
#' @param col_sd Name (in quotes) of the column with within-lab precision.
#' @param LoB Limit of blank (LoB).
#' @param N Total number of measurement values per sample per reagent lot (number of measurements
#' per sample times number of sample). If all reagent lots are pooled, use the total number
#' of measurement values per sample (sum across all lots).
#' @param model Choice of precision model. Default is to choose the model with the lowest AIC.
#' One can also select to use a linear, quadratic, or Sadler model.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns the limit of detection (LoD) value as calculated with the precision
#' profile approach.
#'
#' @examples
#' reagent_lot <- c(rep(1, 6), rep(2, 6))
#' sample <- rep(c(A,B,C,D,E,F), times = 2)
#' avg <- c(.69, 1.42, 2.65, 4.08, 6.08, 10.36, .78, 1.73, 2.89, 3.82, 6.33, 10.92)
#' sd_wl <- c(.39, .39, .46, .55, .64, 1.12, .29, .54, .55, .63, .82, 1.38)
#'
#' sample_df <- data.frame(reagent_lot, sample, avg, sd_wl)
#'
#' LoD(sample_df, col_lot = "reagent_lot", col_sample = "sample", col_avg = "avg",
#' col_sd = "sd_wl", LoB = 0.51, N = 80*6)
#'
#' @export

LoD_precision_profile <- function(df, col_lot, col_sample, col_avg, col_sd, LoB, N,
                                  model = c("lowest AIC", "linear", "quadratic", "Sadler"),
                                  beta = 0.05, always_sep_lots = FALSE){

  model <- match.arg(model)

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    col_lot <- "lot_number"
    df[[col_lot]] <- 1
  }

  # rename columns in df
  names(df)[names(df) == col_sd] <- "sd_wl"
  names(df)[names(df) == col_avg] <- "avg"

  # percentile
  pct <- 1 - beta

  # find number of reagent lots
  n_lots <- unique(df[[col_lot]]) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df[[col_lot]] <- 1
    n_lots <- 1
  }

  # warning about always_sep_lots
  if(always_sep_lots & n_lots > 3){
    message("Warning: Since there are at least four reagent lots in the data provided, CLSI guidelines recommend combining all reagent lots. Consider setting `always_sep_lots` = FALSE.")
  }

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df[[col_lot]])


  #----- precision models -----#

  mod_names <- c("linear", "quadratic", "sadler")

  # linear
  lin_mod <- lapply(lots_list, function(x) lm(sd_wl ~ avg, data = x))

  # quadratic
  quad_mod <- lapply(lots_list, function(x) lm(sd_wl ~ avg + I(avg^2), data = x))

  # Sadler
  # sadler_mod <- lapply(lots_list, function(x)
  #   nls(sd_wl ~ I((a + b*avg)), start = list(a = 0, b = 1), data = sample_df))
  sadler_mod <- lin_mod
  # fix this later!

  all_mods <- list(lin_mod, quad_mod, sadler_mod) |> setNames(mod_names)

  # find model with lowest AIC for each reagent lot
  min_AIC <- vector()
  for(l in 1:n_lots){
    min_AIC[l] <- c(lapply(lin_mod, AIC)[[l]],
                    lapply(quad_mod, AIC)[[l]],
                    lapply(sadler_mod, AIC)[[l]]) |> which.min()
  }

  # which model is necessary? (most complex of all reagent lots)
  mod_nb <- max(min_AIC)
  best_mod <- mod_names[mod_nb]

  # if user has selected a model
  if(model != "lowest AIC"){
    final_mod <- model
    # warning if this model does not have the lowest AIC
    if(best_mod != final_mod){
      message(paste0("Warning: The ", best_mod, " model has a lower AIC than the selected model."))
    }
  # if user wants model with lowest AIC
  }else{
    final_mod <- best_mod
  }

  #----- plot precision profiles -----#

  # set x and y limits for plot
  x_lims <- c(min(sample_df[[col_avg]] - 0.1), max(sample_df[[col_avg]]) + 0.1)
  y_lims <- c(min(sample_df[[col_sd]] - 0.1), max(sample_df[[col_sd]]) + 0.1)

  # plot profile for first reagent lot
  lot_1 <- df[df[[col_lot]] == 1,]
  plot(lot_1[[col_avg]], lot_1[[col_sd]],
       main = "Check model fit!", xlab = "Measurand", ylab = "Within-lab Precision", col = 2,
       xlim = x_lims, ylim = y_lims, pch = 16)
  # and precision model
  new_vals <- seq(from = x_lims[1], to = x_lims[2], length = 20)
  new_preds <- predict(all_mods[[final_mod]][[1]], new_vals |> as.data.frame() |> setNames("avg"))
  lines(new_vals, new_preds, col = 2, lty = 2)

  # plot remaining profiles (if any)
  if(n_lots > 1){
    # loop through remaining reagent lots
    for(l in 2:n_lots){
      lot_l <- df[df[[col_lot]] == l,]
      points(lot_l[[col_avg]], lot_l[[col_sd]], col = l+1, pch = 16)
      # and precision model
      new_preds <- predict(all_mods[[final_mod]][[l]], new_vals |> as.data.frame() |> setNames("avg"))
      lines(new_vals, new_preds, col = l+1, lty = 2)
    }
  }

  #----- calculate LoD -----#

  # # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df[[col_lot]] == l,]

    # critical value
    N_tot <- N # total number of measurements per reagent lot
    K <- unique(lot_l[[col_sample]]) |> length() # number of samples
    cp <- qnorm(pct) / (1 - (1/(4*(N_tot-K)) ))

    LoD_vals[l] <- LoB + cp*predict(all_mods[[final_mod]][[l]],
                                    LoB |> as.data.frame() |> setNames("avg"))
  }

  # LoD val is correct for lot 1 but incorrect for lot 2!!
  # weird, Nick got the same as me...



  # look at all lots together


  # return(LoD_vals)
}
