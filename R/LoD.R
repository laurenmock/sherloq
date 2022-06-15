
#' Limit of Detection (LoD) Calculation
#'
#' @param df A data frame with low-level samples.
#' @param col_lot Name (in quotes) of the column with reagent lot number.
#' @param col_value Name (in quotes) of the column with measurements.
#' @param LoB Limit of blank (LoB). This value should be the same for all reagent lots.
#' @param approach Method of LoD calculation.
#' @param beta Type II error. Default is 0.05.
#' @param always_sep_lots If FALSE, reagent lots are evaluated according to CLSI guidelines
#' (all lots evaluated separately if 2 or 3 lots, and all lots evaluated together if >3 lots).
#' If TRUE, all reagent lots are evaluated separately regardless of the number of lots.
#' Default is FALSE.
#'
#' @return Returns the limit of detection (LoD) value(s).
#'
#' @example
#' reagent_lot <- c(rep(1, 12*5), rep(2, 12*5))
#' day <- rep(rep(c(1, 2, 3), each = 4, times = 10))
#' sample <- rep(c(1, 2, 3, 4, 5), each = 12, times = 2)
#' replicate <- rep(c(1, 2, 3, 4), times = 3*5*2)
#' pg.ml <- c(21.0, 22.8, 28.2, 25.9, 26.4, 28.3, 20.7, 21.9,
#'             24.7, 22.5, 28.5, 29.2, 13.3, 12.6, 18.2, 14.7,
#'             17.8, 14.0, 14.1, 12.5, 11.3, 12.2, 16.2, 13.9,
#'             12.8, 12.9, 17.4, 16.0, 15.9, 14.1, 11.3, 9.4,
#'             10.6, 13.6, 17.6, 14.9, 17.3, 19.2, 21.5, 22.2,
#'             24.1, 25.8, 16.0, 16.4, 24.9, 23.8, 22.1, 26.1,
#'             19.2, 22.7, 28.3, 26.2, 25.1, 30.3, 23.4, 19.2,
#'             26.3, 23.1, 27.5, 30.1, 22.0, 22.5, 21.8, 22.1,
#'             20.3, 21.0, 25.3, 26.0, 27.2, 25.1, 25.3, 25.3,
#'             15.6, 21.2, 14.8, 14.9, 16.0, 15.8, 21.6, 22.8,
#'             15.3, 18.7, 18.3, 19.5, 13.0, 15.9, 9.0, 7.0,
#'             13.4, 8.5, 16.3, 18.1, 12.4, 11.1, 11.3, 10.1,
#'             18.8, 17.6, 14.1, 14.9, 19.2, 15.8, 19.8, 21.4,
#'             18.0, 18.0, 19.6, 23.1, 32.9, 30.4, 29.4, 27.6,
#'             27.7, 30.6, 31.4, 30.4, 32.5, 28.9, 29.8, 35.1))
#' lowlvl_df <- data.frame(reagent_lot, day, sample, replicate, pg.ml)
#'
#' LoD(lowlvl_df, col_lot = "reagent_lot", col_value = "pg.ml", LoB = 7.1, approach = "classical")
#'
#' @export

LoD <- function(df, col_lot, col_value, LoB,
                approach = c("classical", "precision profile", "probit"),
                beta = 0.05, always_sep_lots = FALSE){

  # percentile
  pct <- 1 - beta

  # find number of reagent lots
  n_lots <- unique(df[[col_lot]]) |> length()

  # if there are more than 3 lots and always_sep_lots = FALSE OR all reagent lot values are the same,
  # reset all reagent lot values to 1
  if((n_lots > 3 & !always_sep_lots) | n_lots == 1){
    df[[col_lot]] <- 1
    n_lots <- 1
  }

  if(always_sep_lots & n_lots > 3){
    message("Warning: Since there are at least four reagent lots in the data provided, CLSI guidelines recommend combining all reagent lots. Consider setting `always_sep_lots` = FALSE.")
  }

  # loop through each reagent lot
  LoD_vals <- list()
  for(l in 1:n_lots){

    # look at each lot separately
    lot_l <- df[df[[col_lot]] == l,]

    # classical (parametric)
    if(approach == "classical"){

      # test to see if normally distributed
      shapiro_pval <- shapiro.test(df[[col_value]])$p.value |> round(4)
      if(shapiro_pval <= 0.05){
        message(paste0("Warning: These values do not appear to be normally distributed (Shapiro-Wilk test p-value = ", shapiro_pval,
                       "). Consider a mathematical transformation or a different approach."))
      }

      # standard deviation for each sample in each reagent lot
      sd_sample_df <- do.call(data.frame, aggregate(pg.ml ~ sample, data = lot_l,
                                                    FUN = function(x) c(sd = sd(x), n = length(x))))
      # pooled standard deviation
      sd_L <- with(sd_sample_df, sqrt( ((pg.ml.n-1) %*% pg.ml.sd^2) / sum(pg.ml.n-1)) )

      # critical value
      L <- nrow(lot_l) # number of results
      J <- unique(lot_l$sample) |> length() # number of samples
      cp <- qnorm(pct) / (1 - (1/(4*(L-J)) ))

      # calculate LoD
      LoD_vals[l] <- LoB + cp*sd_L
      names(LoD_vals)[l] <- paste0("LoD_lot_", l)
    }

    # precision profile
    if(approach == "precision profile"){

      # mean and standard deviation (SD_WL) for each sample in each reagent lot
      mean_SDwl_df <- do.call(data.frame, aggregate(pg.ml ~ sample, data = lot_l,
                                                    FUN = function(x) c(xbar = mean(x), sd = sd(x))))
      # not just sd! SD_WL (within-laboratory precision)

      # plot precision profiles
      plot(mean_SDwl_df$pg.ml.xbar, mean_SDwl_df$pg.ml.sd,
           xlab = "Measurand", ylab = "Within-lab Precision")

      # fit first-order polynomial
      pp_model1 <- lm(pg.ml.sd ~ pg.ml.xbar, data = mean_SDwl_df)
      AIC1 <- AIC(pp_model1)
      # fit second-order polynomial
      pp_model2 <- lm(pg.ml.sd ~ poly(pg.ml.xbar, 2, raw = TRUE), data = mean_SDwl_df)
      AIC2 <- AIC(pp_model2)

      # critical value
      N_tot <- nrow(lot_l)
      K <- unique(lot_l$sample) |> length() # number of samples
      cp <- qnorm(pct) / (1 - (1/(4*(N_tot-K)) ))

      # select model with lower AIC
      if(AIC1 <= AIC2){
        init_LoD_vals[l] <- LoB + cp*(c(1, .51) %*% summary(pp_model1)$coef[,1])
      }else{
        init_LoD_vals[l] <- LoB + cp*(c(1, .51, .51^2) %*% summary(pp_model2)$coef[,1])
      }


    }


  }

  # look at all lots from precision profile
  if(approach == "precision profile"){

  }

  # if multiple LoD values, find the max
  if(length(LoD_vals) > 1){
    LoD_vals[n_lots + 1] <- unlist(LoD_vals) |> max()
    names(LoD_vals)[n_lots + 1] <- "LoD_max"
  }

  return(LoD_vals)
}
