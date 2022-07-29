
# note: I haven't written tests to check model selection/starting coefficients, because
# this needs to be updated

#--- data for tests ---#

LoD_PP_df <- data.frame(Reagent = c(rep(1, 6), rep(2, 6)),
                        Sample = rep(LETTERS[1:6], 2),
                        Mean = c(.69, 1.42, 2.65, 4.08, 6.08, 10.36,
                                 .78, 1.73, 2.89, 3.82, 6.33, 10.92),
                        SD_within_lab = c(.39, .39, .46, .55, .64, 1.12,
                                          .29, .54, .55, .63, .82, 1.38))

# column with 4 reagent lots
LoD_PP_df$Reagent4 <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))

# measurement column with a missing value
LoD_PP_df$Mean_NA <- c(NA, LoD_PP_df$Mean[2:nrow(LoD_PP_df)])


#--- tests ---#

suppressMessages(test_that("LoD (precision profile) values are correct", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$lot_1,
               1.167,
               tolerance = 0.001)
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$lot_2,
               1.166,
               tolerance = 0.001)
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$reported,
               1.167,
               tolerance = 0.001)
}))


suppressMessages(test_that("can provide n_samples instead of col_sample", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     n_samples = 6)$LoD_values$reported,
               1.167,
               tolerance = 0.001)
}))


suppressMessages(test_that("can set col_lot = NULL to pool all lots", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = NULL,
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$reported,
               1.174,
               tolerance = 0.001)
}))


suppressMessages(test_that("properly pools lots when there are >3 lots", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent4",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$reported,
               1.174,
               tolerance = 0.001)
}))


suppressMessages(test_that("can change beta", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     beta = 0.01,
                                     LoB = 0.51,
                                     n_measures = 80*6)$LoD_values$reported,
               1.5,
               tolerance = 0.001)
}))


suppressMessages(test_that("warning if >3 lots and always_sep_lots = TRUE", {
  expect_warning(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent4",
                                       col_sample = "Sample",
                                       col_avg = "Mean",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6,
                                       always_sep_lots = TRUE))
}))


test_that("error if column names aren't in df", {
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "lot",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6))
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "lot",
                                     col_sample = "sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6))
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "lot",
                                     col_sample = "Sample",
                                     col_avg = "avg",
                                     col_sd = "sd_wl",
                                     LoB = 0.51,
                                     n_measures = 80*6))
})


suppressMessages(test_that("warning if col_sample and n_samples are both provided", {
  expect_warning(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent",
                                       col_sample = "Sample",
                                       col_avg = "Mean",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6,
                                       n_samples = 6))
}))


suppressMessages(test_that("error if neither col_sample nor n_samples are provided", {
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent",
                                       col_avg = "Mean",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6))
}))


suppressMessages(test_that("warning if missing data in relevant columns", {
  expect_warning(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent",
                                       col_sample = "Sample",
                                       col_avg = "Mean_NA",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6))
}))


suppressMessages(test_that("can select specific model", {
  expect_equal(LoD_precision_profile(df = LoD_PP_df |> subset(Reagent == 2),
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     model = "linear")$LoD_values$reported,
               1.13,
               tolerance = 0.001)
  expect_equal(LoD_precision_profile(df = LoD_PP_df |> subset(Reagent == 1),
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     model = "quadratic")$LoD_values$reported,
               1.167,
               tolerance = 0.001)
  expect_equal(LoD_precision_profile(df = LoD_PP_df |> subset(Sample != "F"),
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     model = "sadler",
                                     sadler_start = c(1,1,1))$LoD_values$reported,
               1.163,
               tolerance = 0.001)
}))


suppressMessages(test_that("error if we select Sadler without starting coefficients", {
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent",
                                       col_sample = "Sample",
                                       col_avg = "Mean",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6,
                                       model = "sadler"))
}))


suppressMessages(test_that("error if sadler_start is the wrong length", {
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     model = "sadler",
                                     sadler_start = c(1,1)))
}))


suppressMessages(test_that("error if sadler_start coefficients don't work", {
  expect_error(LoD_precision_profile(df = LoD_PP_df,
                                     col_lot = "Reagent",
                                     col_sample = "Sample",
                                     col_avg = "Mean",
                                     col_sd = "SD_within_lab",
                                     LoB = 0.51,
                                     n_measures = 80*6,
                                     model = "sadler",
                                     sadler_start = c(1,1,1)))
}))


suppressMessages(test_that("warning if the selected model doesn't have the lowest AIC", {
  expect_warning(LoD_precision_profile(df = LoD_PP_df,
                                       col_lot = "Reagent",
                                       col_sample = "Sample",
                                       col_avg = "Mean",
                                       col_sd = "SD_within_lab",
                                       LoB = 0.51,
                                       n_measures = 80*6,
                                       model = "linear"))
}))
