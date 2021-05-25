
message("   Critical Values")


test_that("Critical values are calculated correctly", {

  ccrit_1 <- fracdist_values(iq = 1, iscon = 0, bb = 0.73, ipc = FALSE, clevel = 0.05)
  ccrit_1_test <- 3.7066

  expect_equal(ccrit_1, ccrit_1_test)

  ccrit_2 <- fracdist_values(iq = 3, iscon = 1, bb = 1.27, ipc = FALSE, clevel = 0.05)
  ccrit_2_test <- 38.7691

  expect_equal(ccrit_2, ccrit_2_test)

  ccrit_3 <- fracdist_values(iq = 12, iscon = 1, bb = 1.27, ipc = FALSE, clevel = 0.05)
  ccrit_3_test <- 407.8485

  expect_equal(ccrit_3, ccrit_3_test)

})


