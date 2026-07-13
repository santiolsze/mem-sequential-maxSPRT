#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))

test_dir("experimentacion/shiny/tests/testthat", reporter = "summary")
test_dir("experimentacion/tests/testthat", reporter = "summary")
