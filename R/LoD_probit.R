
#' Limit of Detection (LoD) Calculation--Probit Approach
#'
#' When should I use the probit approach to calculate LoD?
#'
#' From CLSI EP17 guidelines:
#' Molecular measurement procedures (eg, for infectious diseases, nucleic acid testing)
#' differ from typical measurement procedures because all blank or negative sample
#' results normally are reported as negative. The false-positive rate is much lower
#' than 5% (typically below 0.5%), and the LoB, typically, is taken to be zero. The LoD
#' is calculated from a probit regression model as the measurand concentration at
#' which, with a predefined probability (usually 95%), measurement results yield a
#' positive classification.
#'
#' CLSI EP17 Requirements:
#' - Two reagent lots
#' - One instrument system
#' - Three days
#' - Three positive samples
#' - 30 individual negative patient samples
#' - Five dilutions per positive sample
#' - 20 replicates per dilution (across all testing days) per positive sample per reagent lot
#' - Two replicates (across all testing days) per negative sample per reagent lot
#'
#' @param df This function takes two possible data frame shapes. The data frame can contain
#' either 1) 0s/1s representing negative/positive results for each row or 2) a column with each
#' unique concentration level, a column with the total number of observed positives at that level,
#' and a column with the number of total calls at that level.
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (default).
#' To split the LoD results by any other variable (e.g. lab), simply include the name
#' of this other variable here and set always_sep_lots = TRUE.
#' @param col_conc Name (in quotes) of the column with the concentration.
#' @param col_01 Name (in quotes) of the column with 0s and 1s for negative and positive results.
#' Should be NULL if col_obs_pos and col_tot are provided.
#' @param col_obs_pos Name (in quotes) of the column with the number of positive results. Should be
#' NULL if col_01 is provided.
#' @param col_tot Name (in quotes) of the column with the total number of results. Should be NULL
#' if col_01 is provided.
#' @param LoB Limit of blank (LoB). If the false positives in each reagent lot do not exceed
#' (100*alpha)%, the LoB can be set to zero. Otherwise, calculate the LoB using the classical
#' approach.
#' @param log10_trans According to CLSI EP17 guidelines, a log 10 transformation of concentration
#' values can lead to a better probit fit. If FALSE, a note will indicate if the model(s) with a log
#' 10 transformation of concentration values have a better fit.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns a list with the fitted probit model(s), a plot with the model(s),
#' the model predictions that were used to create the plot, the LoD values for each reagent
#' lot (if evaluated separately), and the reported overall LoD.
#'
#' @examples
#' # CLSI EP17 Appendix C
#' Concentration <- rep(c(0, .006, .014, .025, .05, .15, .3, .5), each = 3)
#' Reagent <- rep(c(1, 2, 3), times = 8)
#' Observed_Positives <- c(0, 0, 0, 11, 12, 22, 15, 22, 31, 23, 28, 27, 29,
#' 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32)
#' Total_Calls <- c(22, 22, 22,30, 30, 34, 30, 30, 34, rep(32, 15))
#'
#' LoD_P_df <- data.frame(Concentration, Reagent, Observed_Positives, Total_Calls)
#'
#' results <- LoD_probit(df = LoD_P_df,
#'                       col_lot = "Reagent",
#'                       col_conc = "Concentration",
#'                       col_obs_pos = "Observed_Positives",
#'                       col_tot = "Total_Calls",
#'                       LoB = 0)
#'
#' @export


LoD_probit <- function(df, col_lot = NULL, col_conc,
                       col_01 = NULL, col_obs_pos = NULL, col_tot = NULL,
                       LoB, log10_trans = FALSE, beta = 0.05, always_sep_lots = FALSE){

  # confirm that column names exist in df
  stopifnot("`col_lot` is not a column in df" = col_lot %in% names(df))
  stopifnot("`col_01` is not a column in df" = col_01 %in% names(df))
  stopifnot("`col_conc` is not a column in df" = col_conc %in% names(df))
  stopifnot("`col_obs_pos` is not a column in df" = col_obs_pos %in% names(df))
  stopifnot("`col_tot` is not a column in df" = col_tot %in% names(df))

  # check for missing data
  relev_cols <- c(col_lot, col_conc, col_01, col_obs_pos, col_tot)
  if(!all(stats::complete.cases(df[,relev_cols]))){
    # remove rows with missing values (and give warning)
    df <- df[stats::complete.cases(df),]
    warning("Ignoring rows with missing values.")
  }

  # if column for reagent lot is NULL, make a column with a vector of 1s (all lot 1)
  if(is.null(col_lot)){
    df$lot <- 1
  }

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_conc] <- "conc"
  names(df)[names(df) == col_01] <- "call"
  names(df)[names(df) == col_obs_pos] <- "obs_pos"
  names(df)[names(df) == col_tot] <- "tot"

  # confirm that user provided col_01 or col_obs_pos and col_tot
  if(is.null(col_01)){
    stopifnot("since col_01 is NULL, col_obs_pos and col_tot must both be provided" =
                !is.null(col_obs_pos) & !is.null(col_tot))
  # if col_01 is provided, then col_obs_pos and col_tot should not be provided
  } else if(!is.null(col_obs_pos) | !is.null(col_tot)){
    stop("col_01 is provided, so col_obs_pos and col_tot should not be provided.
              Please see documentation (?LoD_probit) for more details.")
  # otherwise confirm that col_01 is all 0s and 1s
  } else {
    stopifnot("col_01 must contain only 0s and 1s" = all(df$call %in% c(0,1)))
  }

  # confirm that columns are numeric
  stopifnot("`col_lot` must be numeric" = is.numeric(df$lot))
  stopifnot("`col_conc` must be numeric" = is.numeric(df$conc))
  if(is.null(col_01)){
    stopifnot("`col_obs_pos` must be numeric" = is.numeric(df$obs_pos))
    stopifnot("`col_tot` must be numeric" = is.numeric(df$tot))
  }

  # if user has provided data that is not yet aggregated -- aggregate it!
  if(!is.null(col_01)){
    df <- do.call(data.frame,
                  stats::aggregate(df$call ~ df$conc + df$lot,
                                   data = df,
                                   FUN = function(x) c(obs_pos = sum(x), n = length(x)))) |>
      stats::setNames(c("conc", "lot", "obs_pos", "tot"))
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

  # find hit rate
  df$hit_rate <- df$obs_pos / df$tot

  # remove rows with concentration of 0
  df <- df[df$conc != 0,]

  # log transform concentration vals
  df$conc_log10 <- df$conc |> log10()

  # make reagent lots separate elements in a list
  lots_list <- split(df, f = df$lot)

  # fit GLMs with 1) concentrations and 2) log10 concentrations

  mod_glm <- lapply(lots_list, function(x)
    suppressWarnings(stats::glm(cbind(obs_pos, tot - obs_pos) ~ conc, data = x,
                         family = stats::binomial(link = "probit"))))

  mod_glm_log10 <- lapply(lots_list, function(x)
    suppressWarnings(stats::glm(cbind(obs_pos, tot - obs_pos) ~ conc_log10, data = x,
                         family = stats::binomial(link = "probit"))))

  # get model AICs
  AIC_glm <- sapply(mod_glm, stats::AIC)
  AIC_glm_log10 <- sapply(mod_glm_log10, stats::AIC)

  # select model based on user input
  if(log10_trans){
    mod <- mod_glm_log10
  }else{
    mod <- mod_glm
  }

  # get model coefficients
  # mod_coeff <- lapply(mod, function(x) summary(x)$coef[,1])
  # names(mod_coeff) <- paste0("lot_", 1:n_lots)

  # check model fit with deviance
  chisq_p_val <- lapply(mod, function(x)
    1 - stats::pchisq(sum(stats::residuals(x, type = "pearson")^2), df = x$df.residual))

  # if log transformed
  if(log10_trans){
    new_vals <- seq(log(0.0001), max(df$conc_log10), by = 0.001) # possible log conc. values
    new_preds <- lapply(mod, function(x) stats::predict(x,
                                                 newdata = data.frame("conc_log10" = new_vals),
                                                 se.fit = TRUE)[1:2] |> as.data.frame())
    pred_df <- do.call(rbind, new_preds)
    pred_df$lot <- rep(1:n_lots, each = length(new_vals))
    pred_df$conc <- 10^(new_vals)

    # if not log transformed
  }else{
    new_vals <- seq(0, max(df$conc), by = 0.0001) # possible conc. values
    new_preds <- lapply(mod, function(x) stats::predict(x,
                                                 newdata = data.frame("conc" = new_vals),
                                                 se.fit = TRUE)[1:2] |> as.data.frame())
    pred_df <- do.call(rbind, new_preds)
    pred_df$lot <- rep(1:n_lots, each = length(new_vals))
    pred_df$conc <- new_vals
  }

  # get inverse link function
  inverse_link <- stats::family(mod[[1]])$linkinv

  # use inverse link to get predicted probabilities and 95% CI
  pred_df$p_hat <- inverse_link(pred_df$fit)
  pred_df$lwr <- inverse_link(pred_df$fit - 1.96*pred_df$se.fit)
  pred_df$upr <- inverse_link(pred_df$fit + 1.96*pred_df$se.fit)


  #---- plot -----#
  # graphics::plot.new()
  graphics::par(mfrow = c(1, n_lots),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)
  m <- matrix(c(1:n_lots, rep(n_lots+1, n_lots)), nrow = 2, ncol = n_lots, byrow = TRUE)
  graphics::layout(mat = m, heights = c(0.9, 0.1))

  LoD_vals <- list()

  # loop through reagent lots
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df$lot == l,]
    pred_l <- pred_df[pred_df$lot == l,]

    # plot observed data
    with(lot_l,
         plot(conc, hit_rate, type = "p", log = "x", pch = 16,
              main = ifelse(n_lots == 1, "", paste0("Reagent Lot ", l)),
              xlab = "Concentration (log scale)",
              ylab = "Hit Rate",
              ylim = c(0,1),
              axes = FALSE))

    # function to hide unnecessary GLM warning
    hide_warning <- function(w){
      if((any(grepl("fitted probabilities numerically 0 or 1 occurred", w)))){
        invokeRestart("muffleWarning")}
    }

    # use tryCatch to identify if the GLM doesn't converge
    # if it doesn't converge, don't plot model fit
    tryCatch(

      expr = {

        # test to see if the selected model converges
        # hide warning about some probabilities being 0 or 1
        if(log10_trans){
          test_mod <- withCallingHandlers(stats::glm(cbind(obs_pos, tot - obs_pos) ~ conc_log10,
                                              data = lot_l,
                                              family = stats::binomial(link = "probit")),
                                          warning = hide_warning)
        }else{
          test_mod <- withCallingHandlers(stats::glm(cbind(obs_pos, tot - obs_pos) ~ conc,
                                              data = lot_l,
                                              family = stats::binomial(link = "probit")),
                                          warning = hide_warning)
        }

        # shading
        graphics::polygon(c(pred_l$conc, rev(pred_l$conc)), c(pred_l$upr, rev(pred_l$lwr)),
                col = "lightgoldenrod1", border = "black")
        # model fitted values
        graphics::lines(pred_l$conc,
              pred_l$p_hat,
              type="l", col="black")
        # line at 1-beta (typically 0.95) % hit rate
        graphics::abline(h = pct, col = 'red', lty = 2, lwd = 1.8)
        # add points again so they're on top
        graphics::points(lot_l$conc, lot_l$hit_rate, pch = 16)
        graphics::axis(side = 1)
        graphics::axis(side = 2)
      },

      warning = function(w){
        warning("Original error message: ", w,
                "The GLM algorithm did not converge for lot ", l,
                ", likely because most hit rates
                are very close to 0 or 1. More measurements are needed in the range of
                concentration values for which the hit rate is between 0 and 1.")
        # text("GLM didn't\n converge", cex = 0.7,  col = "red",
        #      x = quantile(pred_df$conc, 0.1), y = 0.1)
      }
    )

    # find first predicted probability greater than pct (typically 95%) and set as LoD
    LoD_vals[l] <- pred_l[pred_l$p_hat >= pct,]$conc[1]
    names(LoD_vals)[l] <- paste0("lot_", l)
  }

  # add legend
  graphics::par(mar = c(0,0,0,0))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  graphics::legend(x = "bottom", inset = 0,
         legend = c("95% Confidence Interval"),
         col = "lightgoldenrod1",
         lwd = 10,
         cex = 1, horiz = TRUE, bty = "n")

  # save plot
  hit_rate_plot <- grDevices::recordPlot()
  graphics::par(mfrow = c(1,1))

  # compare model AIC (log transformed vs. not log transformed) for all lots
  # if any reagent lots have better model fit with log transform, use log transform
  log10_better <- any((AIC_glm - AIC_glm_log10) > 0)

  # message for user
  if(!log10_trans & log10_better){
    warning("The probit models have a better fit on at least one reagent lot if a log transformation
    of the concentration values is performed. Consider setting `log10_trans` = TRUE.")
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

  # warning about GOF
  if(any(chisq_p_val < 0.05)){
    if(length(chisq_p_val) == 1){
      warning("Pearson chi-square goodness-of-fit tests indicate that the probit model fit
            may be insufficient.")
    } else {
      bad_fit <- which(chisq_p_val < 0.05)
      warning("Pearson chi-square goodness-of-fit tests indicate that the probit model fit
            may be insufficient for the following reagent lot(s): ",
              bad_fit)
    }
  }

  # add names to make output easier to read
  names(mod) <- paste0("lot_", 1:n_lots)

  output <- list(mod, hit_rate_plot, pred_df, LoD_vals)
  names(output) <- c("probit_model", "hit_rate_plot", "model_predictions", "LoD_values")

  return(output)
}
