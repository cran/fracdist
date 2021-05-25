
message("   Low Fractional Integration Order")


test_that("P-values are calculated correctly", {

  pval_1 <- fracdist_values(iq = 1, iscon = 0, bb = 0.43, stat = 3.84)
  pval_1_test <- 1 - pchisq(q = 3.84, df = 1^2)

  expect_equal(pval_1, pval_1_test)

  pval_2 <- fracdist_values(iq = 3, iscon = 1, bb = 0.27, stat = 32.84)
  pval_2_test <- 1 - pchisq(q = 32.84, df = 3^2)

  expect_equal(pval_2, pval_2_test)

  pval_3 <- fracdist_values(iq = 12, iscon = 1, bb = 0.27, stat = 412.84)
  pval_3_test <- 1 - pchisq(q = 412.84, df = 12^2)

  expect_equal(pval_3, pval_3_test)

})


test_that("Critical values are calculated correctly", {

  ccrit_1 <- fracdist_values(iq = 1, iscon = 0, bb = 0.43, ipc = FALSE, clevel = 0.05)
  ccrit_1_test <- qchisq(p = 1 - 0.05, df = 1^2)

  expect_equal(ccrit_1, ccrit_1_test)

  ccrit_2 <- fracdist_values(iq = 3, iscon = 1, bb = 0.27, ipc = FALSE, clevel = 0.05)
  ccrit_2_test <- qchisq(p = 1 - 0.05, df = 3^2)

  expect_equal(ccrit_2, ccrit_2_test)

  ccrit_3 <- fracdist_values(iq = 12, iscon = 1, bb = 0.27, ipc = FALSE, clevel = 0.05)
  ccrit_3_test <- qchisq(p = 1 - 0.05, df = 12^2)

  expect_equal(ccrit_3, ccrit_3_test)

})


