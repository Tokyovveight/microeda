if (!requireNamespace("testthat", quietly = TRUE)) {
  message("Skipping tests because testthat is not installed.")
  quit(save = "no", status = 0)
}

library(testthat)
library(microeda)

test_check("microeda")
