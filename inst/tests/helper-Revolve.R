library(testthat)
library(Revolve)

is_greater_than <- function(value) {
  function(actual)
    expectation(actual > value, paste("is not greater than", value))
}

is_less_than <- function(value) {
  function(actual)
    expectation(actual < value, paste("is not less than", value))
}
