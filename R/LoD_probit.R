
#' Limit of Detection (LoD) Calculation--Probit Approach
#'
#' @param df A data frame with...
#' @param col_lot Name (in quotes) of the column with reagent lot number. Can be NULL (no quotes).
#' @param col_conc Name (in quotes) of the column with the concentration.
#' @param col_obs_pos Name (in quotes) of the column with the number of positive results.
#' @param col_tot Name (in quotes) of the column with the total number of results.
#' @param LoB Limit of blank (LoB). Can be set to zero if it can be confirmed that LoB = 0, or
#' calculated using the classical approch.
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

  # confirm that col_conc, col_obs_pos, and col_sd are numeric and positive
  stopifnot("`col_conc` must be numeric" =
              is.numeric(df[[col_conc]])) & all(df[[col_conc]] >= 0)
  stopifnot("`col_obs_pos` must be numeric" =
              is.numeric(df[[col_obs_pos]])) & all(df[[col_obs_pos]] >= 0)
  stopifnot("`col_tot` must be numeric" =
              (is.numeric(df[[col_tot]]) & all(df[[col_tot]] > 0)))

  # percentile
  pct <- 1 - beta

  # find number of reagent lots
  n_lots <- unique(df[[col_lot]]) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE, reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots)){
    df[[col_lot]] <- 1
    n_lots <- 1
  }

  # rename columns in df
  names(df)[names(df) == col_conc] <- "concentration"
  names(df)[names(df) == col_obs_pos] <- "obs_positives"
  names(df)[names(df) == col_tot] <- "tot_calls"


  # find hit rate
  df[["hit_rate"]] <- df[[col_obs_pos]] / df[[col_tot]]

  # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df[[col_lot]] == l,]

    # fit glm
    mod_glm <- glm(cbind(obs_positives, tot_calls - obs_positives) ~ concentration, data = lot_l,
                   family = binomial("probit"))

    inverse_link <- family(mod_glm)$linkinv
    pred_df <- as.data.frame(predict(mod_glm,
                                          newdata = data.frame(concentration = seq(0, .5, length = 100)),
                                          se.fit = TRUE)[1:2])
    pred_df$concentration <- seq(0, .5, length = 100)
    pred_df$p_hat <- inverse_link(pred_df$fit)
    pred_df$lwr <- inverse_link(pred_df$fit - 1.96*pred_df$se.fit)
    pred_df$upr <- inverse_link(pred_df$fit + 1.96*pred_df$se.fit)

    with(lot_l,
         plot(concentration, hit_rate, type = "p",
              xlim = c(0, .5)))
    lines(pred_df$concentration,
          pred_df$p_hat,
          type="l", col="black")
    lines(pred_df$concentration,
          pred_df$upr,
          type="l", col="grey50")
    lines(pred_df$concentration,
          pred_df$lwr,
          type="l", col="grey50")
    abline(h = .95, col = 'red', lty = 2)

    LoD_vals[l] <- pred_df[pred_df$p_hat > pct,]$concentration[1]
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
  return(LoD_vals)
}
