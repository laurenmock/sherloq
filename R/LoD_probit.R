
#' Limit of Detection (LoD) Calculation--Probit Approach
#'
#' @param df A data frame with...
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_conc Name (in quotes) of the column with the concentration.
#' @param col_obs_pos Name (in quotes) of the column with the number of positive results.
#' @param col_tot Name (in quotes) of the column with the total number of results.
#' @param LoB Limit of blank (LoB). If the false positives in each reagent lot do not exceed
#' (100*alpha)%, the LoB can be set to zero. Otherwise, calculate the LoB using the classical approach.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @export

LoD_probit <- function(df, col_lot, col_conc, col_obs_pos, col_tot,
                       LoB, beta = 0.05, always_sep_lots = FALSE){

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
  stopifnot("`col_conc` is not a column in df" = col_conc %in% names(df))
  stopifnot("`col_obs_pos` is not a column in df" = col_obs_pos %in% names(df))
  stopifnot("`col_tot` is not a column in df" = col_tot %in% names(df))

  # rename columns in df
  names(df)[names(df) == col_lot] <- "lot"
  names(df)[names(df) == col_conc] <- "conc"
  names(df)[names(df) == col_obs_pos] <- "obs_pos"
  names(df)[names(df) == col_tot] <- "tot"

  # confirm that all columns are numeric and positive
  stopifnot("`col_lot` must be numeric" =
              (is.numeric(df$lot) & all(df$lot >= 0)))
  stopifnot("`col_conc` must be numeric" =
              is.numeric(df$conc)) & all(df$conc >= 0)
  stopifnot("`col_obs_pos` must be numeric" =
              is.numeric(df$obs_pos)) & all(df$obs_pos >= 0)
  stopifnot("`col_tot` must be numeric" =
              (is.numeric(df$tot) & all(df$tot >= 0)))

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

  # loop through each reagent lot
  LoD_vals <- list()
  chisq_p_val <- vector()
  chisq_stat <- vector()
  plot.new()
  par(mfrow = c(1, n_lots),
      mar = c(3, 3, 2, 1),
      mgp = c(2, 0.5, 0),
      tck = -.01)
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df$lot == l,]

    # fit glm
    mod_glm <- glm(cbind(obs_pos, tot - obs_pos) ~ conc, data = lot_l,
                   family = binomial(link = "probit"))

    # check model fit with deviance
    # chisq_p_val[l] <- 1 - (pchisq(deviance(mod_glm), mod_glm$df.residual))
    chisq_stat <- sum(residuals(mod_glm, type = "pearson")^2)
    chisq_p_val[l] <- 1 - pchisq(chisq_stat, df = mod_glm$df.residual)

    # if fit is bad, consider log transformation of concentrations (and be sure to remove rows with
    # concentration of 0 first)






    # use probit inverse to get predicted probabilities and CIs for a range of concentrations
    inverse_link <- family(mod_glm)$linkinv
    pred_df <- as.data.frame(predict(mod_glm,
                                     newdata = data.frame(conc = seq(0, .5, length = 100)),
                                     se.fit = TRUE)[1:2])
    pred_df$conc <- seq(0, .5, length = 100)
    pred_df$p_hat <- inverse_link(pred_df$fit)
    pred_df$lwr <- inverse_link(pred_df$fit - 1.96*pred_df$se.fit)
    pred_df$upr <- inverse_link(pred_df$fit + 1.96*pred_df$se.fit)

    # plot curves
    with(lot_l,
         plot(conc, hit_rate, type = "p", log = "x",
              main = paste0("Reagent Lot ", l),
              xlab = "Concentration",
              ylab = "Hit Rate"))
    lines(pred_df$conc,
          pred_df$p_hat,
          type="l", col="black")
    lines(pred_df$conc,
          pred_df$upr,
          type="l", col="grey50")
    lines(pred_df$conc,
          pred_df$lwr,
          type="l", col="grey50")
    abline(h = .95, col = 'red', lty = 2)

    # find first probability over 95% and set as LoD
    LoD_vals[l] <- pred_df[pred_df$p_hat > pct,]$conc[1]
    names(LoD_vals)[l] <- paste0("LoD_lot_", l)

  }
  hit_rate_plot <- recordPlot()
  par(mfrow = c(1,1))

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

  output <- list(LoD_vals, hit_rate_plot, chisq_p_val)
  names(output) <- c("LoD_values", "hit_rate_plot", "GOF")

  return(output)
}
