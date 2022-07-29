
#--- data for tests ---#

# example data
LoB_df <- data.frame(Reagent = c(rep(1, 12*5), rep(2, 12*5)),
                     Day = rep(rep(c(1, 2, 3), each = 4, times = 10)),
                     Sample = rep(LETTERS[1:5], each = 12, times = 2),
                     Replicate = rep(c(1, 2, 3, 4), times = 3*5*2),
                     Measurement = c(2.6, -.8, 5.5, 6.0, 4.5, .6, -2.3, 3.4,
                                     5.9, 7.6, 4.1, -1.4, 1.0, 2.9, 4.9, 8.0,
                                     6.9, 5.0, 3.4, 1.2, 6.5, 5.6, -2.2, 2.3,
                                     -4.4, -3.4, 7.0, 6.9, 4.3, 3.2, -1.4, 4.2,
                                     5.9, 7.6, 3.8, 5.8, 1.5, -1.9, 5.1, 5.7,
                                     4.1, 4.5, -.6, .5, 5.4, 7.6, 4.4, 6.6,
                                     1.2, -.7, 6.1, 5.1, 4.8, 3.3, -2.8, -1.4,
                                     8.7, 3.6, 5.1, 3.5, 4.6, 4.1, 1.6, 3.7,
                                     2.2, .7, 4.6, 2.6, 1.1, -4.4, .9, .7,
                                     9.2, 8.3, 4.8, 5.4, 4.8, 6.3, 5.4, 9.6,
                                     7.7, 3.1, 6.1, 10.0, 6.1, 3.2, 3.9, 1.4,
                                     3.1, 4.1, 1.0, 3.4, .1, .4, 2.9, -1.6,
                                     4.0, 11.5, 4.5, 3.6, 4.4, 6.8, 7.1, 4.2,
                                     3.7, 3.7, 5.3, 4.5, 4.0, 6.2, -.2, 2.3,
                                     1.6, 2.6, 6.4, 5.7, 4.2, 3.7, 1.4, 1.5))

# column with 4 reagent lots
LoB_df$Reagent4 <- c(rep(1,30), rep(2,30), rep(3,30), rep(4,30))

# measurement column with a missing value
LoB_df$Measurement_NA <- c(NA, LoB_df$Measurement[2:nrow(LoB_df)])

#--- tests ---#

test_that("LoB (non-parametric) values are correct", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement")$LoB_values$lot_1,
               7.6,
               tolerance = 0.001)
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement")$LoB_values$lot_2,
               9.4,
               tolerance = 0.001)
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement")$LoB_values$reported,
               9.4,
               tolerance = 0.001)
})


test_that("LoB (parametric) values are correct", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement",
                   parametric = TRUE)$LoB_values$lot_1,
               8.718,
               tolerance = 0.001)
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement",
                   parametric = TRUE)$LoB_values$lot_2,
               8.599,
               tolerance = 0.001)
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement",
                   parametric = TRUE)$LoB_values$reported,
               8.718,
               tolerance = 0.001)
})


test_that("can set col_lot = NULL to pool all lots", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = NULL,
                   col_sample = "Sample",
                   col_value = "Measurement")$LoB_values$reported,
               8.15,
               tolerance = 0.001)
})


test_that("properly pools lots when there are >3 lots", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent4",
                   col_sample = "Sample",
                   col_value = "Measurement")$LoB_values$reported,
               8.15,
               tolerance = 0.001)
})


test_that("can change alpha", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement",
                   alpha = 0.01)$LoB_values$reported,
               11.35,
               tolerance = 0.001)
})


test_that("can switch plot to histogram", {
  expect_equal(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "Measurement",
                   plot = "histogram")$LoB_values$reported,
               9.4,
               tolerance = 0.001)
})


test_that("warning if >3 lots and always_sep_lots = TRUE", {
  expect_warning(LoB(df = LoB_df,
                     col_lot = "Reagent4",
                     col_sample = "Sample",
                     col_value = "Measurement",
                     always_sep_lots = TRUE))
})


test_that("error if column names aren't in df", {
  expect_error(LoB(df = LoB_df,
                   col_lot = "lot",
                   col_sample = "Sample",
                   col_value = "Measurement"))
  expect_error(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "samp",
                   col_value = "Measurement"))
  expect_error(LoB(df = LoB_df,
                   col_lot = "Reagent",
                   col_sample = "Sample",
                   col_value = "measure"))
})


test_that("warning if missing data in relevant columns", {
  expect_warning(LoB(df = LoB_df,
                     col_lot = "Reagent",
                     col_sample = "Sample",
                     col_value = "Measurement_NA"))
})


