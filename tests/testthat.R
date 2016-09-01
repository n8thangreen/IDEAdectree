library(testthat)
library(IDEAdectree)

# test_check("IDEAdectree")
# test_file("./tests/testthat/tester.R")
test_dir("./tests/testthat/", reporter = "summary")
