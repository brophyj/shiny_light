# Top-level testthat runner: finds R/ helpers relative to the repo root
# and runs the tests in tests/testthat.
library(testthat)

# Resolve the repo root whether run from the repo root or the tests dir.
root <- local({
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
  if (is.null(here) || !nzchar(here)) here <- getwd()
  for (cand in c(here, dirname(here), getwd())) {
    if (file.exists(file.path(cand, "R", "priors.R"))) return(cand)
  }
  getwd()
})

Sys.setenv(HSF_ROOT = root)
test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
