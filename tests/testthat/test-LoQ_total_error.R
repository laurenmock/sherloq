
#--- data for tests ---#

LoQ_TE_df <- data.frame(Reagent = c(rep(1, 9*5), rep(2, 9*5)),
                        Day = rep(rep(c(1, 2, 3), each = 3, times = 10)),
                        Sample = rep(c(1, 2, 3, 4, 5), each = 9, times = 2),
                        Reference_Value = rep(c(38.2, 47.1, 44.7, 36.5, 42.8), each = 9, times = 2),
                        Replicate = rep(c(1, 2, 3), times = 3*5*2),
                        Measurement = c(36.7, 37.9, 38.3, 36.8, 33.5, 39.2, 41.3, 37.9, 34.9,
                                        49.9, 50, 48.1, 47.8, 43.9, 45.6, 45.4, 51.5, 45.8,
                                        46.1, 43.1, 39.4, 47.3, 45.8, 44.8, 44.6, 47.3, 38.9,
                                        33.3, 34.2, 34.5, 43.1, 34, 37.1, 35.3, 32.4, 36, 42.9,
                                        41.8, 43.8, 46.3, 43.3, 46, 42.6, 41.4, 42.8, 38.5, 41,
                                        43.2, 36.8, 42.1, 35.8, 36.8, 44.1, 39.5, 45.8, 47.8,
                                        46.6, 46.9, 51.3, 50.5, 44.3, 47.5, 52.4, 46.7, 43.6,
                                        42.4, 46.5, 47.9, 42.7, 42.1, 43.4, 44.7, 35.5, 40, 34,
                                        32.9, 33.1, 38.6, 36.2, 41.4, 33, 42, 44.1, 43.2, 46.6,
                                        45.5, 43.5, 41.4, 48.2, 45.7))

LoQ_TE_bad_df <- data.frame(Reagent = c(rep(1, 12*5), rep(2, 12*5)),
                            Day = rep(rep(c(1, 2, 3), each = 4, times = 10)),
                            Sample = rep(c(1, 2, 3, 4, 5), each = 12, times = 2),
                            Reference_Value = rep(c(26.1, 16.9, 13.1, 20.4, 27.8), each = 12, times = 2),
                            Replicate = rep(c(1, 2, 3, 4), times = 3*5*2),
                            Measurement = c(21.0, 22.8, 28.2, 25.9, 26.4, 28.3, 20.7, 21.9,
                                            24.7, 22.5, 28.5, 29.2, 13.3, 12.6, 18.2, 14.7,
                                            17.8, 14.0, 14.1, 12.5, 11.3, 12.2, 16.2, 13.9,
                                            12.8, 12.9, 17.4, 16.0, 15.9, 14.1, 11.3, 9.4,
                                            10.6, 13.6, 17.6, 14.9, 17.3, 19.2, 21.5, 22.2,
                                            24.1, 25.8, 16.0, 16.4, 24.9, 23.8, 22.1, 26.1,
                                            19.2, 22.7, 28.3, 26.2, 25.1, 30.3, 23.4, 19.2,
                                            26.3, 23.1, 27.5, 30.1, 22.0, 22.5, 21.8, 22.1,
                                            20.3, 21.0, 25.3, 26.0, 27.2, 25.1, 25.3, 25.3,
                                            15.6, 21.2, 14.8, 14.9, 16.0, 15.8, 21.6, 22.8,
                                            15.3, 18.7, 18.3, 19.5, 13.0, 15.9, 9.0, 7.0,
                                            13.4, 8.5, 16.3, 18.1, 12.4, 11.1, 11.3, 10.1,
                                            18.8, 17.6, 14.1, 14.9, 19.2, 15.8, 19.8, 21.4,
                                            18.0, 18.0, 19.6, 23.1, 32.9, 30.4, 29.4, 27.6,
                                            27.7, 30.6, 31.4, 30.4, 32.5, 28.9, 29.8, 35.1))


# column with 5 reagent lots
LoQ_TE_df$Reagent5 <- c(rep(1, 18), rep(2, 18), rep(3, 18), rep(4, 18), rep(5, 18))

# measurement column with a missing value
LoQ_TE_df$Measurement_NA <- c(NA, LoQ_TE_df$Measurement[2:nrow(LoQ_TE_df)])


#--- tests ---#

test_that("LoQ (total error) values are correct", {
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6)$LoQ_values$lot_1,
               35.54,
               tolerance = 0.001)
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6)$LoQ_values$lot_2,
               36.08,
               tolerance = 0.001)
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6)$LoQ_values$reported,
               36.08,
               tolerance = 0.001)
})


test_that("can plot linear model fit instead of connecting points", {
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6,
                               plot_lm = TRUE)$LoQ_values$reported,
               36.08,
               tolerance = 0.001)
})


test_that("error if all points are above the accuracy goal", {
  expect_error(LoQ_total_error(df = LoQ_TE_bad_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6))
})


test_that("can set col_lot = NULL to pool all lots", {
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = NULL,
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6)$LoQ_values$reported,
               35.81,
               tolerance = 0.001)
})


test_that("properly pools lots when there are >3 lots", {
  expect_equal(LoQ_total_error(df = LoQ_TE_df,
                               col_lot = "Reagent5",
                               col_sample = "Sample",
                               col_ref = "Reference_Value",
                               col_value = "Measurement",
                               accuracy_goal = 21.6)$LoQ_values$reported,
               35.81,
               tolerance = 0.001)
})


test_that("warning if >3 lots and always_sep_lots = TRUE", {
  expect_warning(LoQ_total_error(df = LoQ_TE_df,
                                 col_lot = "Reagent5",
                                 col_sample = "Sample",
                                 col_ref = "Reference_Value",
                                 col_value = "Measurement",
                                 accuracy_goal = 21.6,
                                 always_sep_lots = TRUE))
})


test_that("error if column names aren't in df", {
  expect_error(expect_equal(LoQ_total_error(df = LoQ_TE_df,
                                            col_lot = "lot",
                                            col_sample = "Sample",
                                            col_ref = "Reference_Value",
                                            col_value = "Measurement",
                                            accuracy_goal = 21.6)))
  expect_error(expect_equal(LoQ_total_error(df = LoQ_TE_df,
                                            col_lot = "Reagent",
                                            col_sample = "sample",
                                            col_ref = "Reference_Value",
                                            col_value = "Measurement",
                                            accuracy_goal = 21.6)))
  expect_error(expect_equal(LoQ_total_error(df = LoQ_TE_df,
                                            col_lot = "Reagent",
                                            col_sample = "Sample",
                                            col_ref = "ref",
                                            col_value = "Measurement",
                                            accuracy_goal = 21.6)))
  expect_error(expect_equal(LoQ_total_error(df = LoQ_TE_df,
                                            col_lot = "Reagent",
                                            col_sample = "Sample",
                                            col_ref = "Reference_Value",
                                            col_value = "measure",
                                            accuracy_goal = 21.6)))
})


test_that("warning if missing data in relevant columns only", {
  expect_warning(LoQ_total_error(df = LoQ_TE_df,
                                 col_lot = "Reagent",
                                 col_sample = "Sample",
                                 col_ref = "Reference_Value",
                                 col_value = "Measurement_NA",
                                 accuracy_goal = 21.6))
})

