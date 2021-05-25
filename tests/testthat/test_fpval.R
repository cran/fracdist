
message("   P-values")


test_that("P-values are calculated correctly", {

  pval_1 <- fracdist_values(iq = 1, iscon = 0, bb = 0.73, stat = 3.84)
  pval_1_test <- 0.0461

  expect_equal(pval_1, pval_1_test)

  pval_2 <- fracdist_values(iq = 3, iscon = 1, bb = 1.27, stat = 32.84)
  pval_2_test <- 0.1882

  expect_equal(pval_2, pval_2_test)

  pval_3 <- fracdist_values(iq = 12, iscon = 1, bb = 1.27, stat = 412.84)
  pval_3_test <- 0.0320

  expect_equal(pval_3, pval_3_test)

})


