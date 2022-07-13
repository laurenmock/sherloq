
#' Limit of Detection (LoD) Calculation--Precision Profile Approach
#'
#' When should I use the precision profile approach to calculate LoD?
#'
#' From CLSI EP17 guidelines:
#' The precision profile approach is useful when the variability of measurement results
#' changes significantly in the region of the assumed LoD. It also may be useful in cases
#' in which the developer wishes to make use of precision data that were acquired over a
#' wider concentration interval than typically used for detection capability studies.
#' Implementation of this approach assumes that variability of measurement results vs
#' the mean measurand concentration can be fitted adequately with a precision profile model.
#'
#' @param df A data frame with with the mean value and within-laboratory precision (obtained
#' from CLSI EP05) for each sample from each reagent lot.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_avg Name (in quotes) of the column with measurements.
#' @param col_sd_wl Name (in quotes) of the column with within-lab precision.
#' @param LoB_val Limit of blank (LoB).
#' @param N Total number of measurement values per sample per reagent lot (number of measurements
#' per sample times number of sample). If all reagent lots are pooled, use the total number
#' of measurement values per sample (sum across all lots).
#' @param model Choice of precision model. Default is to choose either the linear or quadratic model
#' based on AIC. One can also select to use a linear, quadratic, or Sadler model, the three models
#' most widely used in clinical literature (according to CLSI guidelines).
#' @param sadler_start Vector of length 3 with starting coefficient values for the Sadler model,
#' with the form (beta0 + beta1*x)^beta2.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns a list with the limit of detection (LoD) value as calculated with the precision
#' profile approach, precision model coefficients, and precision profile plot.
#'
#' @examples
#' # CLSI EP17 Appendix B
#' reagent_lot <- c(rep(1, 6), rep(2, 6))
#' sample <- rep(c("A","B","C","D","E","F"), times = 2)
#' avg <- c(.69, 1.42, 2.65, 4.08, 6.08, 10.36, .78, 1.73, 2.89, 3.82, 6.33, 10.92)
#' sd_wl <- c(.39, .39, .46, .55, .64, 1.12, .29, .54, .55, .63, .82, 1.38)
#'
#' LoD_PP_df <- data.frame(reagent_lot, sample, avg, sd_wl)
#'
#' LoD_precision_profile(df = LoD_PP_df, col_lot = "reagent_lot", col_sample = "sample",
#' col_avg = "avg", col_sd_wl = "sd_wl", LoB = 0.51, N = 80*6)
#'
#' @export

LoD_precision_profile <- function(df, col_lot, col_sample, col_avg, col_sd_wl, LoB_val, N,
                                  model = c("lowest AIC", "linear", "quadratic", "sadler"),
                                  sadler_start = NULL, beta = 0.05, always_sep_lots = FALSE){

  model <- match.arg(model)

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

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_avg` is not a column in df" = col_avg %in% names(df))
  stopifnot("`col_sd_wl` is not a column in df" = col_sd_wl %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_sample] <- "sample"
  names(df)[names(df) == col_avg] <- "avg"
  names(df)[names(df) == col_sd_wl] <- "sd_wl"

  # confirm that columns are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_avg` must be numeric" = is.numeric(df$avg))
  stopifnot("`col_sd_wl` must be numeric" = is.numeric(df$sd_wl))

  # if model is sadler, check sadler_start for starting values
  if(model == "sadler" & !is.null(sadler_start)){
    stopifnot("`sadler_start` must be NULL or a vector with three values (one for each
    coefficient in the Sadler model)" = length(sadler_start) == 3)
  }

  # percentile
  pct <- 1 - beta

  # find number of reagent lots
  n_lots <- unique(df$lot) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df$lot <- 1
    n_lots <- 1
  }

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df$lot)


  #----- precision models -----#

  # linear
  lin_mod <- lapply(lots_list, function(x) stats::lm(sd_wl ~ avg, data = x))

  # quadratic
  quad_mod <- lapply(lots_list, function(x) stats::lm(sd_wl ~ avg + I(avg^2), data = x))

  mod_names <- c("linear", "quadratic")
  all_mods <- list(lin_mod, quad_mod) |> stats::setNames(mod_names)

  # Sadler
  # Sadler is finicky! only run if Sadler is selected by user
  if(model == "sadler"){

    tryCatch(
      expr = {

        # if user specified starting values
        if(!is.null(sadler_start)){

          sadler_mod <- lapply(lots_list, function(x)
            stats::nls(sd_wl ~ I((c0 + c1*avg)^c2), data = x,
                       start = list(c0 = sadler_start[1],
                                    c1 = sadler_start[2],
                                    c2 = sadler_start[3])))

        # if user didn't specify starting values
        }else{

          # Sadler model form:
          # sd_wl = (c0 + c1*avg)^c2
          # come up with an initial guess for c2
          # take the log of both sides
          # log(sd_wl) = c2*log(c0 + c1*avg)

          # sd_wl^(1/c2) = c0 + c1*avg
          # now fit a linear model
          mods_lm <- lapply(lots_list, function(x) stats::lm(x$sd_wl^(1/1) ~ x$avg))

          # get intercepts and slopes
          ints <- sapply(mods_lm, function(x) summary(x)$coeff[1,1] |> exp())
          slopes <- sapply(mods_lm, function(x) summary(x)$coeff[2,1])

          # fit NLS models
          sadler_mod <- lapply(1:n_lots, function(l)
            stats::nls(sd_wl ~ I((c0 + c1*avg)^c2), data = lots_list[[l]],
                       start = list(c0 = ints[l], c1 = slopes[l], c2 = 0.1)))


        }

        # add Sadler to list of models to compare with AIC
        mod_names[3] <- "sadler"
        all_mods[[3]] <- sadler_mod
        names(all_mods) <- mod_names

      },

      error = function(e){
        stop(paste0("Unable to fit the Sadler model. Try changing the starting values or
        select a different model.
        Original error message: ", e))
      }
    )
  }


  # find model with best (lowest) AIC for each reagent lot
  min_AIC <- vector()
  for(l in 1:n_lots){
    min_AIC[l] <- lapply(all_mods, function(x) x[[l]] |> stats::AIC()) |> which.min()
  }

  # which model is necessary? (most complex of all reagent lots)
  mod_nb <- max(min_AIC)
  best_mod <- mod_names[mod_nb]


  #----- messages about AIC/model fit -----#

  # note about which models are being compared
  if(model == "sadler"){
    message("Note: Comparing linear, quadratic, and Sadler model fits.")
  }else{
    message("Note: Comparing linear and quadratic fits. Select `model` = 'sadler' in order to
            compare linear, quadratic, and Sadler model fits.")
  }


  # if user has selected a specific model (not lowest AIC)
  if(model != "lowest AIC"){
    final_mod <- model

    # if the selected model also has the lowest AIC
    if(best_mod == final_mod){
      message(paste0("The ", best_mod, " model, selected by the user, has the best model fit,
                     as determined by AIC."))
    # warning if the selected model does not have the lowest AIC
    }else{
      warning(paste0("The ", best_mod, " model may have a better model fit
                   than the ", model, " model, as determined by AIC."))
    }

  # if user wants model with lowest AIC
  }else{
    final_mod <- best_mod
    message(paste0("Selecting the ", best_mod, " model, which has a better model fit than the ",
                   mod_names[mod_names != best_mod], " model, as determined by AIC."))
  }

  # get model coefficients to report
  mod_coeff <- lapply(all_mods[[final_mod]], function(x) summary(x)$coef[,1])
  names(mod_coeff) <- paste0("lot_", 1:n_lots)


  #----- plot precision profiles -----#

  # set x and y limits for plot
  x_lims <- c(min(min(df$avg), LoB_val), max(df$avg))
  y_lims <- c(min(df$sd_wl), max(df$sd_wl))

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # plot profiles
  #opar <- par()
  graphics::plot.new()
  graphics::par(mfrow = c(1,1),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)
  plot(1, type = "n",
       main = "Check model fit", xlab = "Measurand", ylab = "Within-lab Precision", col = pal[1],
       xlim = x_lims, ylim = y_lims, pch = 16)

  # new x values to plug into model
  new_vals <- seq(from = x_lims[1], to = x_lims[2], length = 20)
  # model predictions
  new_preds <- lapply(1:n_lots, function(x) stats::predict(all_mods[[final_mod]][[x]], new_vals |>
                                                      as.data.frame() |> stats::setNames("avg")))

  # save model predictions in a data frame
  pred_df <- unlist(new_preds) |> as.data.frame() |> stats::setNames("sd_pred")
  pred_df$lot <- rep(1:n_lots, each = length(new_vals))
  pred_df$avg <- rep(new_vals, times = n_lots)
  pred_df <- pred_df[,c(2,1,3)]
  pred_df <- pred_df[,c("lot", "avg", "sd_pred")]


  # draw observed data and model predictions
  temp <- lapply(1:n_lots, function(x) graphics::points(lots_list[[x]]$avg, lots_list[[x]]$sd_wl,
                                              col = pal[x], pch = 16))
  temp <- lapply(1:n_lots, function(x) graphics::lines(new_vals, new_preds[[x]],
                                             col = pal[x], lty = 2))

  # add legend
  if(n_lots > 1){
    labs <- paste0("Lot ", 1:n_lots)
    graphics::legend("topleft", pch = 16, lty = 2, legend = labs, col = pal, cex = 0.7, bty = "n")
  }

  # save plot
  mod_plot <- grDevices::recordPlot()
  #par(opar)


  #----- calculate LoD -----#

  # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df$lot == l,]

    # critical value
    N_tot <- N # total number of measurements per reagent lot
    K <- unique(lot_l$sample) |> length() # number of samples
    cp <- stats::qnorm(pct) / (1 - (1/(4*(N_tot-K)) ))

    # predict within-lab precision across various measurement concentrations
    trial_mc <- seq(from = LoB_val, to = stats::median(df$avg), by = .01)
    trial_sd <- stats::predict(all_mods[[final_mod]][[l]],
                          newdata = trial_mc |> as.data.frame() |> stats::setNames("avg"))

    # caculate LoD and bias for each trial
    trial_LoD <- LoB_val + cp*trial_sd
    bias <- trial_mc - trial_LoD

    # find measurand with bias closest to 0 and use those LoD values
    LoD_vals[l] <- trial_mc[which.min(abs(bias))]
    names(LoD_vals)[l] <- paste0("lot_", l)
  }

  # reporting LoD
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

  output <- list(mod_coeff, pred_df, mod_plot, LoD_vals)
  names(output) <- c("model_coeff", "model_predictions", "mod_plot", "LoD_values")

  return(output)
}
