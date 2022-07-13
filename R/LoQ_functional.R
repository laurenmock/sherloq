
#' Limit of Quantitation (LoQ) Calculation--Functional Sensitivity Approach
#'
#' @param df A data frame with...
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_sample Name (in quotes) of the column with sample number.
#' @param col_avg Name (in quotes) of the column with measurements.
#' @param col_sd_wl Name (in quotes) of the column with within-lab precision.
#' @param target_cv Desired coefficient of variation (sd/mean)*100.
#' @param model Select the power function (default) or power function with intercept to model
#' coefficient of variation (CV) as a function of measurement value.
#' @param coeff_start Optional vector with starting coefficient values for the selected model.
#' Should have length 2 for the power function, or length 3 for the power function with intercept.
#' Default is NULL, meaning start values will be estimated automatically.
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
#' # CLSI EP17 Appendix D1
#' reagent_lot <- rep(c(1, 2), each = 9)
#' sample <- rep(1:9, times = 2)
#' avg <- c(.04, .053, .08, .111, .137, .164, .190, .214, .245,
#' .041, .047, .077, .106, .136, .159, .182, .205, .234)
#' sd_wl <- c(.016, .016, .016, .017, .014, .012, .011, .016, .013,
#' .018, .014, .012, .019, .016, .015, .015, .016, .014)
#'
#' LoQ_F_df <- data.frame(reagent_lot, sample, avg, sd_wl)
#'
#' LoQ_functional(df = LoQ_F_df, col_lot = "reagent_lot", col_sample = "sample",
#' col_avg = "avg", col_sd_wl = "sd_wl", target_cv = 10, model = "power")
#'
#' @export

LoQ_functional <- function(df, col_lot, col_sample, col_avg, col_sd_wl, target_cv,
                           model = c("power", "power w intercept"),
                           coeff_start = NULL, always_sep_lots = FALSE){

  model <- match.arg(model)

  # if coeff_start estimates were provided, check model with length of initial model coefficients
  if(!is.null(coeff_start)){

    stopifnot("`coeff_start` must be a numeric vector" = is.numeric(coeff_start))

    if(model == "power"){
      stopifnot("The power function requires 2 initial coefficient estimates.
                Set coeff_start = NULL or provide 2 estimates." =
                  length(coeff_start) == 2)
    }else if(model == "power w intercept"){
      stopifnot("The power function with intercept requires 3 initial coefficient estimates.
                Set coeff_start = NULL or provide 3 estimates." =
                  length(coeff_start) == 3)
    }
  }

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

  # confirm that col_avg and col_sd_wl are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_sample` must be numeric" = is.numeric(df$sample))
  stopifnot("`col_avg` must be numeric" = is.numeric(df$avg))
  stopifnot("`col_sd_wl` must be numeric" = (is.numeric(df$sd_wl) & all(df$sd_wl >= 0)))

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



  #----- fit model(s) -----#


  # might need a tryCatch here (telling user to set start_coeff as null if their values don't work)



  # if user didn't provide starting values for the power model:
  if(model == "power" & is.null(coeff_start)){

    # power function:
    # cv = c0*avg^c1
    # take the log of both sides of the power function to get:
    # log(cv) = log(c0) + c1*log(avg)
    # fit a linear model with this new model
    # the intercept of the linear model is log(c0) and the slope is c1
    # so in the original power function model, c0 = exp(int) and c1 = slope
    mods_lm <- lapply(lots_list, function(x) stats::lm(log(x$cv) ~ log(x$avg)))

    # get intercepts and slopes
    ints <- sapply(mods_lm, function(x) summary(x)$coeff[1,1] |> exp())
    slopes <- sapply(mods_lm, function(x) summary(x)$coeff[2,1])

    # fit NLS models
    mods <- lapply(1:n_lots, function(l) stats::nls(cv ~ I(c0*avg^c1), data = lots_list[[l]],
                                                    start = list(c0 = ints[l], c1 = slopes[l])))

  # otherwise just fit power model with coeff_start
  }else if(model == "power"){

    mods <- lapply(lots_list, function(x) nls(cv ~ I(c0*avg^c1), data = x,
                                              start = list(c0 = coeff_start[1],
                                                           c1 = coeff_start[2])))

  # if user didn't provide starting values for the power model with intercept:
  }else if(is.null(coeff_start) & model == "power w intercept"){

    # power function with intercept:
    # cv = c0 + c1*avg^c2
    # move intercept to the left and then take the log of both sides:
    # log(cv - c0) = log(c1) + c2*log(avg)
    # get an initial estimate of c0: half of the CV at the highest measurement value
    cv_at_max <- sapply(lots_list, function(x) x$cv[which.max(x$avg)]/2)
    # fit a linear model with this new model
    # the intercept of the linear model is log(c1) and the slope is c2
    # so in the original power function model, c0 = our initial guess, c1 = exp(int) and c2 = slope
    mods_lm <- lapply(1:n_lots, function(x) stats::lm(log(lots_list[[x]]$cv - cv_at_max[x]) ~
                                                        log(lots_list[[x]]$avg)))

    # get intercepts and slopes
    ints <- sapply(mods_lm, function(x) summary(x)$coeff[1,1] |> exp())
    slopes <- sapply(mods_lm, function(x) summary(x)$coeff[2,1])

    # fit NLS models
    mods <- lapply(1:n_lots, function(l) stats::nls(cv ~ I(c0 + c1*avg^c2), data = lots_list[[l]],
                                                    start = list(c0 = cv_at_min[l],
                                                                 c1 = ints[l],
                                                                 c2 = slopes[l])))

  # otherwise just fit power model with intercept with coeff_start
  }else if(model == "power w intercept"){

    mods <- lapply(lots_list, function(x) stats::nls(cv ~ I(c0 + c1*avg^c1), data = x,
                                                     start = list(c0 = coeff_start[1],
                                                                  c1 = coeff_start[2],
                                                                  c2 = coeff_start[3])))
  }


  # get model coefficients
  mods_coef <- lapply(mods, function(x) summary(x)$coef[,1])

  # calculate LoQ for each lot
  if(model == "power"){
    LoQ_vals <- lapply(mods_coef, function(x) (target_cv/x[1]) ^ (1/x[2]))
  }else{
    LoQ_vals <- lapply(mods_coef, function(x) ((target_cv - x[1])/x[2]) ^ (1/x[3]))
  }


  # naming
  names(mods_coef) <- paste0("lot_", 1:n_lots)
  names(LoQ_vals) <- paste0("lot_", 1:n_lots)

  # remove label c0 (label for each element within each element of the list)
  for(l in 1:n_lots){
    names(LoQ_vals[l][[1]]) <- NULL
  }


  #----- plot -----#

  # color blind safe palette
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  graphics::plot.new()
  graphics::par(mfrow = c(1,1),
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
    new_preds[[l]] <- stats::predict(mods[[l]], new_vals |> as.data.frame() |> stats::setNames("avg"))
    graphics::lines(new_vals, new_preds[[l]], col = pal[l], lty = 2)
  }

  # add legend
  if(n_lots > 1){
    labs <- paste0("Lot ", 1:n_lots)
    graphics::legend("topright", pch = 16, lty = 2, legend = labs, col = pal, cex = 0.7, bty = "n")
  }

  # save plot
  mod_plot <- grDevices::recordPlot()

  # reporting LoQ
  if(always_sep_lots & length(LoQ_vals) > 3){
    warning("Since there are at least four reagent lots in the data provided, CLSI guidelines
            recommend combining all reagent lots. Set `always_sep_lots` = FALSE to obtain a single,
            reportable estimate of LoQ.")

  # if only one LoQ value, report as LoQ_reported (not LoQ_lot_1)
  }else if(length(LoQ_vals) == 1){
    names(LoQ_vals)[1] <- "reported"

  # otherwise find max LoQ to report
  }else{
    LoQ_vals[n_lots + 1] <- unlist(LoQ_vals) |> max()
    names(LoQ_vals)[n_lots + 1] <- "reported"
  }

  output <- list(mods_coef, mod_plot, LoQ_vals)
  names(output) <- c("model_coeff", "model_plot", "LoQ_values")

  return(output)
}

