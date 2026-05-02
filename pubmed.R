# pubmed.R
# Thin caller that delegates to R/pubmed.R::fetch_cv_rcts().
#
# The original monolithic query logic has been moved into the reusable
# `fetch_cv_rcts()` helper. This script preserves the "run-as-a-script" entry
# point so that existing workflows (e.g. `Rscript pubmed.R`) still work.

source("R/pubmed.R")

# Defaults mirror the systematic CV-RCT pipeline (PR 2):
#   - 6 target journals
#   - cardiovascular MeSH / keyword scope
#   - RCT publication type
# Adjust `start`/`end` and/or `journals` below as needed.

fetch_cv_rcts(
  journals = c("N Engl J Med", "Lancet", "BMJ", "JAMA",
               "JAMA Cardiol", "Circulation"),
  start    = "2020/01/01",
  end      = "3000/12/31",
  out_csv  = "data/pubmed_citations.csv"
)
