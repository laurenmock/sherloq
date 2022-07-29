
#--- data for tests ---#

# df with 2 lots
LoQ_F_df <- data.frame(Reagent = rep(c(1, 2), each = 9),
                       Sample = rep(1:9, times = 2),
                       Mean = c(.04, .053, .08, .111, .137, .164, .190, .214, .245,
                                .041, .047, .077, .106, .136, .159, .182, .205, .234),
                       SD_within_lab = c(.016, .016, .016, .017, .014, .012, .011, .016, .013,
                                         .018, .014, .012, .019, .016, .015, .015, .016, .014))

# column with 4 reagent lots
LoQ_F_df$Reagent4 <- c(rep(c(1,2), each = 5),
                       rep(c(3,4), each = 4))

# measurement column with a missing value
LoQ_F_df$SD_within_lab_NA <- c(NA, LoQ_F_df$SD_within_lab[2:nrow(LoQ_F_df)])


#--- tests ---#

test_that("LoQ (functional) values are correct", {
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10)$LoQ_values$lot_1,
               0.1458,
               tolerance = 0.001)
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10)$LoQ_values$lot_2,
               0.1468,
               tolerance = 0.001)
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10)$LoQ_values$reported,
               0.1468,
               tolerance = 0.001)
})


test_that("can set col_lot = NULL to pool all lots", {
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = NULL,
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10)$LoQ_values$reported,
               0.1463,
               tolerance = 0.001)
})


test_that("properly pools lots when there are >3 lots", {
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent4",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10)$LoQ_values$reported,
               0.1463,
               tolerance = 0.001)
})


test_that("warning if >3 lots and always_sep_lots = TRUE", {
  expect_warning(LoQ_functional(df = LoQ_F_df,
                                col_lot = "Reagent4",
                                col_avg = "Mean",
                                col_sd_wl = "SD_within_lab",
                                target_cv = 10,
                                always_sep_lots = TRUE))
})


test_that("error if column names aren't in df", {
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "lot",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10))
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "avg",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10))
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "sd",
                              target_cv = 10))
})


test_that("warning if missing data in relevant columns", {
  expect_warning(LoQ_functional(df = LoQ_F_df,
                                col_lot = "Reagent",
                                col_avg = "Mean",
                                col_sd_wl = "SD_within_lab_NA",
                                target_cv = 10))
})


test_that("can select power function with intercept", {
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10,
                              model = "power w intercept")$LoQ_values$reported,
               0.1464,
               tolerance = 0.001)
})


test_that("can input starting coefficients", {
  expect_equal(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10,
                              coeff_start = c(1.2, -1))$LoQ_values$reported,
               0.1468,
               tolerance = 0.001)
})


test_that("error if starting coefficient vector is the wrong length for the model", {
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10,
                              coeff_start = c(1,2,3)))
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10,
                              model = "power w intercept",
                              coeff_start = c(1,2)))
})


test_that("error if the starting coefficients don't work", {
  expect_error(LoQ_functional(df = LoQ_F_df,
                              col_lot = "Reagent",
                              col_avg = "Mean",
                              col_sd_wl = "SD_within_lab",
                              target_cv = 10,
                              coeff_start = c(10,10)))
})


