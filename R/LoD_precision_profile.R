
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
#' CLSI EP17 Requirements:
#' - Two reagent lots
#' - One instrument system
#' - Five days
#' - Five samples
#' - Five replicates per sample (for each reagent lot, day, and instrument system
#' combination)
#' - 40 replicates per sample (across all days and instrument systems) per reagent lot
#'
#' @param df A data frame with with mean measurand concentrations and within-laboratory
#' precision (obtained from CLSI EP05) for each sample from each reagent lot.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (default).
#' To split the LoB results by any other variable (e.g. lab), simply include the name
#' of this other variable here and set always_sep_lots = TRUE.
#' @param col_sample Name (in quotes) of any column with unique values that correspond to each
#' sample. Can be NULL (default) if n_samples is provided. Does not need to be numeric.
#' @param col_avg Name (in quotes) of the column with mean measurements.
#' @param col_sd_wl Name (in quotes) of the column with within-lab precision.
#' @param LoB Limit of blank (LoB). Can be calculated using LoB function on data with
#' blank samples.
#' @param n_measures If there are 1) two or three reagent lots or 2) always_sep_lots = TRUE,
#' meaning that the lots will be evaluated separately, n_measures is the total number of
#' measurement values per reagent lot. If there are 1) four or more reagent lots or 2)
#' col_lot = NULL, meaning that all lots will be pooled, n_measures is the total number of
#' measurement values across all reagent lots.
#' @param n_samples Unique number of samples. Only needs to be provided if col_sample = NULL.
#' @param model Choice of precision model. Default is to choose either the linear or quadratic
#' model based on AIC. You can also select to use a linear, quadratic, or Sadler model, the
#' three models most widely used in clinical literature (according to CLSI guidelines). If
#' you select the Sadler model, you must also provide starting coefficients with sadler_start.
#' @param sadler_start Vector of length 3 with starting coefficient values for the Sadler model,
#' with the form (beta0 + beta1*x)^beta2.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns a list with the fitted precision profile model, the precision profile
#' model plot, model predictions that were used to create the plot, the LoD values for each
#' reagent lot (if evaluated separately), and the reported overall LoD.
#'
#' @examples
#' # CLSI EP17 Appendix B
#' Reagent <- c(rep(1, 6), rep(2, 6))
#' Sample <- rep(c("A","B","C","D","E","F"), times = 2)
#' Mean <- c(.69, 1.42, 2.65, 4.08, 6.08, 10.36, .78, 1.73, 2.89, 3.82, 6.33, 10.92)
#' SD_within_lab <- c(.39, .39, .46, .55, .64, 1.12, .29, .54, .55, .63, .82, 1.38)
#'
#' LoD_PP_df <- data.frame(Reagent, Sample, Mean, SD_within_lab)
#'
#' results <- LoD_precision_profile(df = LoD_PP_df,
#'                                  col_lot = "Reagent",
#'                                  col_sample = "Sample",
#'                                  col_avg = "Mean",
#'                                  col_sd = "SD_within_lab",
#'                                  LoB = 0.51,
#'                                  n_measures = 80*6)
#'
#' @export

LoD_precision_profile <- function(df, col_lot = NULL, col_sample = NULL, col_avg, col_sd_wl,
                                      LoB, n_measures, n_samples = NULL,
                                      model = c("lowestAIC", "linear", "quadratic", "sadler"),
                                      sadler_start = NULL, beta = 0.05, always_sep_lots = FALSE){

  model <- match.arg(model)

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_sample` is not a column in df" = col_sample %in% names(df))
  stopifnot("`col_avg` is not a column in df" = col_avg %in% names(df))
  stopifnot("`col_sd_wl` is not a column in df" = col_sd_wl %in% names(df))

  # check for missing data
  relev_cols <- c(col_lot, col_sample, col_avg, col_sd_wl)
  if(!all(stats::complete.cases(df[,relev_cols]))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # warning if providing col_sample and n_samples
  if(!is.null(col_sample) & !is.null(n_samples)){
    warning("Both `col_sample` and `n_samples` have been provided. Only col_sample will be used." )

  # error if neither col_sample nor n_samples is provided
  } else if(is.null(col_sample) & is.null(n_samples)){
    stop("You must provide either col_sample or n_samples.")
  }

  # make a new column for sample
  if(!is.null(col_sample)){
    df$sample <- df[[col_sample]]
  }

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_avg] <- "Mean"
  names(df)[names(df) == col_sd_wl] <- "SD_within_lab"

  # confirm that columns are numeric
  tryCatch(
    expr = {
      df$lot <- as.numeric(df$lot)
      df$Mean <- as.numeric(df$Mean)
      df$SD_within_lab <- as.numeric(df$SD_within_lab)
    },
    warning = function(w){
      stop("`col_lot`, `col_avg`, and `col_sd_wl` must be numeric")
    }
  )

  # if model is sadler, check sadler_start for starting values
  if(model == "sadler" & !is.null(sadler_start)){
    stopifnot("`sadler_start` must be NULL or a vector with three values (one for each
    coefficient in the Sadler model)" = length(sadler_start) == 3)
  } else if(model == "sadler" & is.null(sadler_start)){
    stop("The Sadler model requires you to specify three starting coefficient values with
    `sadler_start`. You must specify these values or select a different model.")
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
  lin_mod <- lapply(lots_list, function(x) stats::lm(SD_within_lab ~ Mean, data = x))

  # quadratic
  quad_mod <- lapply(lots_list, function(x) stats::lm(SD_within_lab ~ Mean + I(Mean^2), data = x))

  mod_names <- c("linear", "quadratic")
  all_mods <- list(lin_mod, quad_mod) |> stats::setNames(mod_names)

  # Sadler
  # Sadler is finicky! only run if Sadler is selected by user
  if(model == "sadler"){

    #if(is.null(sadler_start)){

      # Sadler model form:
      # cv = (c0 + c1*Mean)^c2
      # write out the exponent
      # cv =
      # take the log of both sides:
      # log(cv) =


    #-------------------------#

    # lot1 <- df |> subset(lot == 1)
    # lot1 <- LoD_PP_df |> subset(Reagent == 1)
    # #lot2 <- df |> subset(lot == 2)
    #
    # model1 <- function(a, data) {
    #   (a[1] + a[2]*data$Mean)^(a[3])
    # }
    #
    # measure_distance <- function(mod, data) {
    #   diff <- data$SD_within_lab - model1(mod, data)
    #   sqrt(mean(diff ^ 2))
    # }
    #
    # best <- optim(c(1, 0, 15), measure_distance, data = lot1)
    # best$par
    #
    # #---------#
    #
    # library(minpack.lm)
    #
    # fit <- nlsLM(SD_within_lab ~ a + b*Mean^c, data = lot1, start = list(a=1, b=1, c=1))
    # summary(fit)
    #
    #
    # #}


    tryCatch(
      expr = {


        sadler_mod <- lapply(lots_list, function(x)
          stats::nls(SD_within_lab ~ I((c0 + c1*Mean)^c2), data = x,
                     start = list(c0 = sadler_start[1],
                                  c1 = sadler_start[2],
                                  c2 = sadler_start[3])))

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
  if(model != "lowestAIC"){
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
  # mod_coeff <- lapply(all_mods[[final_mod]], function(x) summary(x)$coef[,1])
  # names(mod_coeff) <- paste0("lot_", 1:n_lots)


  #----- plot precision profiles -----#

  # set x and y limits for plot
  x_lims <- c(min(min(df$Mean), LoB), max(df$Mean))
  y_lims <- c(min(df$SD_within_lab), max(df$SD_within_lab))

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
       main = "Check model fit", xlab = "Mean measurement", ylab = "Within-lab Precision", col = pal[1],
       xlim = x_lims, ylim = y_lims, pch = 16)

  # new x values to plug into model
  new_vals <- seq(from = x_lims[1], to = x_lims[2], length = 20)
  # model predictions
  new_preds <- lapply(1:n_lots, function(x) stats::predict(all_mods[[final_mod]][[x]], new_vals |>
                                                             as.data.frame() |> stats::setNames("Mean")))

  # save model predictions in a data frame
  pred_df <- unlist(new_preds) |> as.data.frame() |> stats::setNames("sd_pred")
  pred_df$lot <- rep(1:n_lots, each = length(new_vals))
  pred_df$Mean <- rep(new_vals, times = n_lots)
  pred_df <- pred_df[,c(2,1,3)]
  pred_df <- pred_df[,c("lot", "Mean", "sd_pred")]


  # draw observed data and model predictions
  temp <- lapply(1:n_lots, function(x) graphics::points(lots_list[[x]]$Mean, lots_list[[x]]$SD_within_lab,
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

    if(!is.null(col_sample)){
      n_samples <- unique(lot_l$sample) |> length() # number of samples
    }

    # critical value
    cp <- stats::qnorm(pct) / (1 - (1/(4*(n_measures - n_samples)) ))

    # predict within-lab precision across various measurement concentrations
    trial_mc <- seq(from = min(0, LoB), to = max(df$Mean), by = .001)
    trial_sd <- stats::predict(all_mods[[final_mod]][[l]],
                               newdata = trial_mc |> as.data.frame() |> stats::setNames("Mean"))

    # caculate LoD and bias for each trial
    trial_LoD <- LoB + cp*trial_sd
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

  # add names to make output easier to read
  names(all_mods[[final_mod]]) <- paste0("lot_", 1:n_lots)

  # nicer names for pred_df
  names(pred_df) <- c("Reagent", "Mean", "SD_within_lab")

  output <- list(all_mods[[final_mod]], pred_df, mod_plot, LoD_vals)
  names(output) <- c("model", "plot_data", "plot", "LoD_values")

  return(output)
}
