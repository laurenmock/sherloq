
#--- data for tests ---#

# df with 2 lots
LoD_C_df <- data.frame(Reagent = c(rep(1, 12*5), rep(2, 12*5)),
                       Day = rep(rep(c(1, 2, 3), each = 4, times = 10)),
                       Sample = rep(LETTERS[1:5], each = 12, times = 2),
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

# column with 4 reagent lots
LoD_C_df$Reagent4 <- c(rep(1,30), rep(2,30), rep(3,30), rep(4,30))

# measurement column with a missing value
LoD_C_df$Measurement_NA <- c(NA, LoD_C_df$Measurement[2:nrow(LoD_C_df)])


#--- tests ---#

test_that("LoD (classical) values are correct", {
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5)$LoD_values$lot_1,
               12.64,
               tolerance = 0.001)
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5)$LoD_values$lot_2,
               11.94,
               tolerance = 0.001)
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5)$LoD_values$reported,
               12.64,
               tolerance = 0.001)
})


test_that("can set col_lot = NULL to pool all lots", {
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = NULL,
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5)$LoD_values$reported,
               13.00,
               tolerance = 0.001)
})


test_that("properly pools lots when there are >3 lots", {
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent4",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5)$LoD_values$reported,
               13.00,
               tolerance = 0.001)
})


test_that("can change beta", {
  expect_equal(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5,
                             beta = 0.01)$LoD_values$reported,
               14.773,
               tolerance = 0.001)
})


test_that("warning if >3 lots and always_sep_lots = TRUE", {
  expect_warning(LoD_classical(df = LoD_C_df,
                               col_lot = "Reagent4",
                               col_sample = "Sample",
                               col_value = "Measurement",
                               LoB = 7.5,
                               always_sep_lots = TRUE))
})


test_that("error if column names aren't in df", {
  expect_error(LoD_classical(df = LoD_C_df,
                             col_lot = "lot",
                             col_sample = "Sample",
                             col_value = "Measurement",
                             LoB = 7.5))
  expect_error(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "sample",
                             col_value = "Measurement",
                             LoB = 7.5))
  expect_error(LoD_classical(df = LoD_C_df,
                             col_lot = "Reagent",
                             col_sample = "Sample",
                             col_value = "measure",
                             LoB = 7.5))
})


test_that("warning if missing data in relevant columns", {
  expect_warning(LoD_classical(df = LoD_C_df,
                               col_lot = "Reagent",
                               col_sample = "Sample",
                               col_value = "Measurement_NA",
                               LoB = 7.5))
})

