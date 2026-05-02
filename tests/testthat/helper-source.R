# Shared helper: source the R/ modules into the test environment.
local({
  root <- Sys.getenv("HSF_ROOT", unset = "")
  if (!nzchar(root) || !dir.exists(root)) {
    root <- tryCatch(dirname(dirname(dirname(sys.frame(1)$ofile))),
                     error = function(e) getwd())
  }
  for (f in c("priors.R", "likelihoods.R", "eb_mixture.R",
              "summaries.R", "plots.R")) {
    src <- file.path(root, "R", f)
    if (file.exists(src)) source(src, local = FALSE)
  }
})
