context("Test our example")

test_that("our function works",{
  set.seed(100)
  sim1 <- sim_data(1000, 5, 20, 400,alpha0 = 25, alpha_1=25, alpha_2=25)
  #Need to think of better tests
  expect_equal(dim(sim1$X), c(1000,20))
  
})