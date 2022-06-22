
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
#' One can also select to use a linear, quadratic, or Sadler model, the three models most widely used
#' in clinical literature (according to CLSI guidelines).
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

  # check for missing data
  if(!all(complete.cases(df))){
    # remove rows with missing values (and give warning)
    df <- df[complete.cases(df),]
    message("Warning: Ignoring rows with missing values.")
  }

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    col_lot <- "lot_number"
    df[[col_lot]] <- 1
  }

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_avg` is not a column in df" = col_avg %in% names(df))
  stopifnot("`col_sd` is not a column in df" = col_sd %in% names(df))

  # confirm that col_avg and col_sd are numeric
  stopifnot("`col_avg` must be numeric" = is.numeric(df[[col_avg]]))
  stopifnot("`col_sd` must be numeric" = (is.numeric(df[[col_sd]]) & all(df[[col_sd]] >= 0)))

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

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df[[col_lot]])


  #----- precision models -----#

  mod_names <- c("linear", "quadratic", "sadler")

  # linear
  lin_mod <- lapply(lots_list, function(x) lm(sd_wl ~ avg, data = x))

  # quadratic
  quad_mod <- lapply(lots_list, function(x) lm(sd_wl ~ avg + I(avg^2), data = x))

  # Sadler

  # define structure of Sadler precision model (CLSI EP17 5.4.1.3, page 20)
  # sim1 <- data.frame(x = iris$Sepal.Length, y = iris$Sepal.Width, z = iris$Petal.Length)
  # model1 <- function(beta, data){
  #   (beta[1] + data$x * beta[2]) ^ (beta[3])}
  # model1(c(7, 1.5, 1), sim1)
  # measure_distance <- function(mod, data) {
  #   diff <- data$y - model1(mod, data)
  #   sqrt(mean(diff ^ 2))}
  # best <- optim(c(0, 0, 0), measure_distance, data = sim1)
  # best$par
#
#   sadler <- function(beta, data){
#     (beta[1] + beta[2] * data$avg) ^ (beta[3])}
#
#   # sadler(c(1,1,1), sample_df)
#
#   # root mean square error
#   rmse_fun <- function(mod, data) {
#     diff <- data$avg - sadler(mod, data)
#     sqrt(mean(diff^2))}
#
#   # sadler_best <- optim(c(0, 0, 0), rmse_fun, data = sample_df)
#   # sadler_best$par
#
#   sadler_best_1 <- optim(c(1, 1, 1), rmse_fun, data = sample_df |> subset(reagent_lot == 1))
#
#   newdat <- as.data.frame(new_vals) |> setNames("avg")
#   coefs <- sadler_best_1$par
#   new_preds_sad <- sadler(beta = coefs, data = newdat)
#
  # sadler_mod <- lapply(lots_list, function(x))
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

    # if this model does have the lowest AIC
    if(best_mod == final_mod){
      message(paste0("The ", best_mod, " model, selected by the user, has the best model fit as
                     determined by AIC."))
    # warning if this model does not have the lowest AIC
    }else{
      message(paste0("Warning: The ", best_mod, " model may have a better model fit
                     than the selected model, as determined by AIC."))
    }

  # if user wants model with lowest AIC
  }else{
    final_mod <- best_mod
    message(paste0("Selecting the ", best_mod, " model, which has the best model fit
                   of the three choices as determined by AIC."))
  }

  #----- plot precision profiles -----#

  # set x and y limits for plot
  x_lims <- c(min(sample_df[[col_avg]] - 0.1), max(sample_df[[col_avg]]) + 0.1)
  y_lims <- c(min(sample_df[[col_sd]] - 0.1), max(sample_df[[col_sd]]) + 0.1)

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # plot profile for first reagent lot
  lot_1 <- df[df[[col_lot]] == 1,]
  plot.new()
  plot(lot_1[[col_avg]], lot_1[[col_sd]],
       main = "Check model fit", xlab = "Measurand", ylab = "Within-lab Precision", col = pal[1],
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
      points(lot_l[[col_avg]], lot_l[[col_sd]], col = pal[l+1], pch = 16)
      # and precision model
      new_preds <- predict(all_mods[[final_mod]][[l]], new_vals |> as.data.frame() |> setNames("avg"))
      lines(new_vals, new_preds, col = pal[l+1], lty = 2)
    }
  }
  mod_plot <- recordPlot()

  #----- calculate LoD -----#

  # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df[[col_lot]] == l,]

    # critical value
    N_tot <- N # total number of measurements per reagent lot
    K <- unique(lot_l[[col_sample]]) |> length() # number of samples
    cp <- qnorm(pct) / (1 - (1/(4*(N_tot-K)) ))

    init_LoD <- LoB + cp*predict(all_mods[[final_mod]][[l]],
                                  newdata = LoB |> as.data.frame() |> setNames("avg"))

    # now calculate what we would get using different measurand concentrations (instead of 0.51, LoB)
    mc <- seq(from = LoB, to = init_LoD*2, by = .01)
    trial_lod <- LoB + cp*predict(all_mods[[final_mod]][[l]],
                             newdata = mc |> as.data.frame() |> setNames("avg"))
    bias <- mc - trial_lod

    # find measurand with bias closest to 0 and use those LoD values
    LoD_vals[l] <- mc[which.min(abs(bias))]
    names(LoD_vals)[l] <- paste0("LoD_lot_", l)
  }

  # warning about always_sep_lots when n_lots > 3
  if(always_sep_lots & length(LoD_vals) > 3){
    message("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a single,
            reportable estimate of LoD.")
    # if only one LoD value, report as LoB_reported (not LoD_lot_1)
  }else if(length(LoD_vals) == 1){
    names(LoD_vals)[1] <- "LoD_reported"
    # otherwise find max LoD to report
  }else{
    LoD_vals[n_lots + 1] <- unlist(LoD_vals) |> max()
    names(LoD_vals)[n_lots + 1] <- "LoD_reported"
  }

  output <- list(LoD_vals, mod_plot)
  names(output) <- c("LoD_values", "precision_model_plot")

  return(output)
}
