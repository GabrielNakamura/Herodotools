test_that("calc_harmonic_mean works on normal numeric vectors", {
  expect_equal(calc_harmonic_mean(c(1, 2, 3)),
               3 / (1/1 + 1/2 + 1/3))
})

test_that("calc_harmonic_mean handles zeros correctly", {
  # With a zero (no ignoring): should collapse to 0
  expect_equal(calc_harmonic_mean(c(1, 2, 0)), 0)
  
  # With a zero (ignore.zero = TRUE): should ignore the zero
  expect_equal(calc_harmonic_mean(c(1, 2, 0), ignore.zero = TRUE),
               2 / (1/1 + 1/2))
})

test_that("calc_harmonic_mean handles NA values correctly", {
  # With NA (na.rm = FALSE): should return NA
  expect_true(is.na(calc_harmonic_mean(c(1, 2, NA))))
  
  # With NA removed
  expect_equal(calc_harmonic_mean(c(1, 2, NA), na.rm = TRUE),
               2 / (1/1 + 1/2))
  
  # With NA and zero, both removed
  expect_equal(calc_harmonic_mean(c(1, 2, 0, NA), na.rm = TRUE, ignore.zero = TRUE),
               2 / (1/1 + 1/2))
})

test_that("calc_harmonic_mean handles edge cases", {
  # Empty vector
  expect_true(is.na(calc_harmonic_mean(numeric(0))))
  
  # Only zeros, ignore.zero = TRUE
  expect_true(is.na(calc_harmonic_mean(c(0, 0, 0), ignore.zero = TRUE)))
  
  # Only NAs, na.rm = TRUE
  expect_true(is.na(calc_harmonic_mean(c(NA, NA), na.rm = TRUE)))
})

