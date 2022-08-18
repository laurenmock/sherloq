
#' Limit of Quantitation (LoQ) Calculation--Functional Sensitivity Approach
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
#' @param df A data frame with mean measurand concentrations and within-laboratory
#' precision (obtained from CLSI EP05) for each sample from each reagent lot.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (default).
#' To split the LoQ results by any other variable (e.g. lab), simply include the name
#' of this other variable here and set always_sep_lots = TRUE.
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
#' @return Returns a list with the fitted power or power with intercept model(s), a plot with
#' the model(s), the LoQ values for each reagent lot (if evaluated separately), and the
#' reported overall LoQ.
#'
#' @examples
#' # CLSI EP17 Appendix D1
#' Reagent <- rep(c(1, 2), each = 9)
#' Sample <- rep(1:9, times = 2)
#' Mean <- c(.04, .053, .08, .111, .137, .164, .190, .214, .245,
#' .041, .047, .077, .106, .136, .159, .182, .205, .234)
#' SD_within_lab <- c(.016, .016, .016, .017, .014, .012, .011, .016, .013,
#' .018, .014, .012, .019, .016, .015, .015, .016, .014)
#'
#' LoQ_F_df <- data.frame(Reagent, Sample, Mean, SD_within_lab)
#'
#' results <- LoQ_functional(df = LoQ_F_df,
#'                           col_lot = "Reagent",
#'                           col_avg = "Mean",
#'                           col_sd_wl = "SD_within_lab",
#'                           target_cv = 10)
#'
#' @export

LoQ_functional <- function(df, col_lot = NULL, col_avg, col_sd_wl, target_cv,
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
  relev_cols <- c(col_lot, col_avg, col_sd_wl)
  if(!all(stats::complete.cases(df[,relev_cols]))){
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
  stopifnot("`col_avg` is not a column in df" = col_avg %in% names(df))
  stopifnot("`col_sd_wl` is not a column in df" = col_sd_wl %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_avg] <- "avg"
  names(df)[names(df) == col_sd_wl] <- "sd_wl"

  # confirm that columns are numeric
  tryCatch(
    expr = {
      df$lot <- as.numeric(df$lot)
      df$avg <- as.numeric(df$avg)
      df$sd_wl <- as.numeric(df$sd_wl)
    },
    warning = function(w){
      stop("`col_lot`, `col_avg`, and `col_sd_wl` must be numeric")
    }
  )

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


  #----- start plot (will run even if model fit fails) -----#

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


  #----- fit model(s) -----#

  tryCatch(

    expr = {

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

        mods <- lapply(lots_list, function(x) stats::nls(cv ~ I(c0*avg^c1), data = x,
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
                                                        start = list(c0 = cv_at_max[l],
                                                                     c1 = ints[l],
                                                                     c2 = slopes[l])))

        # otherwise just fit power model with intercept with coeff_start
      }else if(model == "power w intercept"){

        mods <- lapply(lots_list, function(x) stats::nls(cv ~ I(c0 + c1*avg^c1), data = x,
                                                         start = list(c0 = coeff_start[1],
                                                                      c1 = coeff_start[2],
                                                                      c2 = coeff_start[3])))
      }

    },

    error = function(e){
      stop(paste0("Original error message: ", e,
                  "Unable to fit the ", model, " model. If you specified starting coefficients
                  with coeff_start, try setting coeff_start = NULL (the default)--the function will
                  estimate starting coefficients for you."))
    }

  )



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


  #----- add model predictions to plot -----#

  new_vals <- seq(0, max(df$avg), length = 100)
  new_preds <- lapply(1:n_lots, function(x)
    stats::predict(mods[[x]], new_vals |> as.data.frame() |> stats::setNames("avg")))
  temp <- lapply(1:n_lots, function(x)
    graphics::lines(new_vals, new_preds[[x]], col = pal[x], lty = 2))

  # add legend
  if(n_lots > 1){
    labs <- paste0("Lot ", 1:n_lots)
    graphics::legend("topright", pch = 16, lty = 2, legend = labs, col = pal, cex = 0.7, bty = "n")
  }

  # save plot
  mod_plot <- grDevices::recordPlot()

  # get data from plot into a table for output
  pred_df <- unlist(new_preds) |> as.data.frame() |> stats::setNames("CV")
  pred_df$Reagent <- rep(1:n_lots, each = length(new_vals))
  pred_df$Measurand <- rep(new_vals, times = n_lots)
  pred_df <- pred_df[,c(2,1,3)]
  pred_df <- pred_df[,c("Reagent", "Measurand", "CV")]

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

  # add names to make output easier to read
  names(mods) <- paste0("lot_", 1:n_lots)

  output <- list(mods, pred_df, mod_plot, LoQ_vals)
  names(output) <- c("model", "plot_data", "plot", "LoQ_values")

  return(output)
}

