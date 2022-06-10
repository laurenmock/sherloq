
#' Limit of Detection (LoD) Calculation
#'
#' @param df data frame for one reagent lot
#' @param col_values Name (in quotes) of the column with measurements.
#' @param beta type II error
#' @param LoB limit of blank
#' @param approach method of LoD calculation
#'
#' @export

# is L really supposed to be across all reagent lots?

LoD <- function(df, col_values, beta=0.05, LoB,
                approach = c("classical", "precision profile", "probit")){

  # # percentile
  # pct <- 1 - alpha
  #
  # # number of reagents
  B <- nrow(df)

  # classical nonparametric
  if(approach == "classical parametric"){

    # maybe?
    # df$pg.ml <- log(df$pg.ml)

    # standard deviation for each sample
    sd_sample_df <- do.call(data.frame,
                             stats::aggregate(pg.ml ~ sample, data = df,
                                       FUN = function(x) c(sd = stats::sd(x), n = length(x))))

    # pooled standard deviation
    sd_L <- with(sd_sample_df, sqrt(((pg.ml.n - 1) %*% pg.ml.sd^2)/sum(pg.ml.n - 1)))

    # critical value
    L <- B # double-check this! says across all reagent lots...
    J <- unique(df$sample) |> length()
    cp <- stats::qnorm(pct) / (1 - (1/(4*(L-J)) ))

    # calculate LoD
    LoD_val <- LoB + cp*sd_L
  }

  if(approach == "precision profile"){
    LoD_val <- 0
  }

}
