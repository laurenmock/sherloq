
#--- data for tests ---#

# df with 3 lots
LoD_P_df <- data.frame(Concentration = rep(c(0, .006, .014, .025, .05, .15, .3, .5), each = 3),
                       Reagent = rep(c(1, 2, 3), times = 8),
                       Observed_Positives = c(0, 0, 0,
                                              11, 12, 22,
                                              15, 22, 31,
                                              23, 28, 27,
                                              29, 32, 32,
                                              32, 32, 32,
                                              32, 32, 32,
                                              32, 32, 32),
                       Total_Calls = c(22, 22, 22,
                                       30, 30, 34,
                                       30, 30, 34,
                                       rep(32, 15)))

# column with 4 reagent lots
LoD_P_df$Reagent4 <- rep(c(1, 2, 3, 4), times = 6)

# measurement column with a missing value
LoD_P_df$Concentration_NA <- c(NA, LoD_P_df$Concentration[2:nrow(LoD_P_df)])

# missing concentrations --> leads to convergence error
LoD_P_df$Concentration_sparse <- ifelse(LoD_P_df$Concentration %in% c(0.006, 0.014),
                                        NA, LoD_P_df$Concentration)

# df that isn't aggregated yet (DEV data)
data(LoDdat_P)

# df that isn't aggregated and has 1s and 2s instead of 0s and 1s
LoDdat_P$SampleCall_bad <- LoDdat_P$SampleCall + 1


#--- tests ---#

test_that("LoD (probit) values are correct", {
  expect_equal(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0)$LoD_values$lot_1,
               0.0561,
               tolerance = 0.001)
  expect_equal(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0)$LoD_values$lot_2,
               0.0299,
               tolerance = 0.001)
  expect_equal(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0)$LoD_values$reported,
               0.0561,
               tolerance = 0.001)
})


test_that("can set col_lot = NULL to pool all lots", {
  expect_warning(LoD_probit(df = LoD_P_df,
                            col_lot = NULL,
                            col_conc = "Concentration",
                            col_obs_pos = "Observed_Positives",
                            col_tot = "Total_Calls",
                            LoB = 0))
})


test_that("properly pools lots when there are >3 lots", {
  expect_equal(suppressWarnings(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent4",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0))$LoD_values$reported,
               0.0425,
               tolerance = 0.001)
})


test_that("can change beta", {
  expect_equal(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0,
                          beta = 0.01)$LoD_values$reported,
               0.0737,
               tolerance = 0.001)
})


test_that("can do a log10 transformation", {
  expect_equal(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0,
                          log10_trans = TRUE)$LoD_values$reported,
               0.0767,
               tolerance = 0.001)
})


test_that("warning if a log10 transformation would be better", {
  w <- capture_warnings((LoD_probit(df = LoD_P_df,
                                    col_lot = "Reagent4",
                                    col_conc = "Concentration",
                                    col_obs_pos = "Observed_Positives",
                                    col_tot = "Total_Calls",
                                    LoB = 0,
                                    always_sep_lots = TRUE)))
  expect_match(w, "log transformation", all = FALSE)
})


test_that("warning if >3 lots and always_sep_lots = TRUE", {
  w <- capture_warnings((LoD_probit(df = LoD_P_df,
                                    col_lot = "Reagent4",
                                    col_conc = "Concentration",
                                    col_obs_pos = "Observed_Positives",
                                    col_tot = "Total_Calls",
                                    LoB = 0,
                                    always_sep_lots = TRUE)))
  expect_match(w, "always_sep_lots", all = FALSE)
})


test_that("warning if the probit fit doesn't pass GOF test", {
  w <- capture_warnings((LoD_probit(df = LoD_P_df,
                                    col_lot = "Reagent4",
                                    col_conc = "Concentration",
                                    col_obs_pos = "Observed_Positives",
                                    col_tot = "Total_Calls",
                                    LoB = 0,
                                    always_sep_lots = TRUE)))
  expect_match(w, "goodness-of-fit", all = FALSE)
})


test_that("warning if at least one model doesn't converge", {
  w <- capture_warnings((LoD_probit(df = LoD_P_df,
                                    col_lot = "Reagent",
                                    col_conc = "Concentration_sparse",
                                    col_obs_pos = "Observed_Positives",
                                    col_tot = "Total_Calls",
                                    LoB = 0)))
  expect_match(w, "converge", all = FALSE)
})


test_that("error if column names aren't in df", {
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "lot",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "conc",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "obs_pos",
                          col_tot = "Total_Calls",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "tot",
                          LoB = 0))
})


test_that("warning if missing data in relevant columns", {
  w <- capture_warnings(LoD_probit(df = LoD_P_df,
                                   col_lot = "Reagent",
                                   col_conc = "Concentration_NA",
                                   col_obs_pos = "Observed_Positives",
                                   col_tot = "Total_Calls",
                                   LoB = 0))
  expect_match(w, "missing", all = FALSE)
})


test_that("error if wrong combination of col_01, col_obs_pos, and col_tot are provided", {
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_tot = "Total_Calls",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_obs_pos = "Observed_Positives",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_01 = "Total_Calls",
                          col_obs_pos = "Observed_Positives",
                          col_tot = "Total_Calls",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          col_01 = "Total_Calls",
                          col_obs_pos = "Observed_Positives",
                          LoB = 0))
  expect_error(LoD_probit(df = LoD_P_df,
                          col_lot = "Reagent",
                          col_conc = "Concentration",
                          LoB = 0))
})


test_that("works on unaggregated data (0s and 1s)", {
  expect_equal(LoD_probit(df = LoDdat_P,
                          col_lot = NULL,
                          col_conc = "MeanObsVAF",
                          col_01 = "SampleCall",
                          LoB = 0)$LoD_values$reported,
               0.0541,
               tolerance = 0.001)
})


test_that("error if col_01 has values other than 0 and 1", {
  expect_error(LoD_probit(df = LoDdat_P,
                          col_lot = NULL,
                          col_conc = "MeanObsVAF",
                          col_01 = "SampleCall_bad",
                          LoB = 0))
})


