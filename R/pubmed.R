# R/pubmed.R
# Refactored PubMed helper for fetching cardiovascular RCTs.
# Depends only on rentrez + xml2 (same dependency choice as the original pubmed.R).

suppressPackageStartupMessages({
  if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
  if (!requireNamespace("xml2",    quietly = TRUE)) install.packages("xml2")
  library(rentrez)
  library(xml2)
})

#' Build a PubMed journal-filter clause
#'
#' Combines `Journal[SB]` / `TA` field tags to restrict a query to the target
#' set of journals. Using both tags is slightly more permissive and tolerant of
#' journal name variants (e.g. NEJM vs. N Engl J Med).
#'
#' @param journals character vector of journal NLM title abbreviations.
#' @return single-string PubMed clause.
.build_journal_clause <- function(journals) {
  stopifnot(length(journals) >= 1)
  pieces <- vapply(journals, function(j) {
    sprintf('("%s"[TA] OR "%s"[Journal])', j, j)
  }, character(1))
  paste0("(", paste(pieces, collapse = " OR "), ")")
}

#' Build the cardiovascular-scope clause (MeSH + publication type).
.build_cv_clause <- function() {
  paste(
    '(',
    '  "cardiovascular diseases"[MeSH Terms]',
    '  OR "heart diseases"[MeSH Terms]',
    '  OR "coronary artery disease"[MeSH Terms]',
    '  OR "myocardial infarction"[MeSH Terms]',
    '  OR "heart failure"[MeSH Terms]',
    '  OR "atrial fibrillation"[MeSH Terms]',
    '  OR "stroke"[MeSH Terms]',
    '  OR "hypertension"[MeSH Terms]',
    '  OR cardiovascular[Title/Abstract]',
    '  OR cardiac[Title/Abstract]',
    ')',
    collapse = "\n"
  )
}

.extract_text <- function(x, xpath) {
  node <- xml_find_first(x, xpath)
  if (inherits(node, "xml_missing")) return(NA_character_)
  xml_text(node)
}

.parse_article <- function(node) {
  pmid     <- .extract_text(node, ".//PMID")
  title    <- .extract_text(node, ".//ArticleTitle")
  journal  <- .extract_text(node, ".//Journal/Title")
  journal_abbrev <- .extract_text(node, ".//Journal/ISOAbbreviation")
  year     <- .extract_text(node, ".//PubDate/Year")
  if (is.na(year)) year <- .extract_text(node, ".//PubDate/MedlineDate")
  abstract <- paste(
    xml_text(xml_find_all(node, ".//Abstract/AbstractText")),
    collapse = " "
  )

  authors <- xml_find_all(node, ".//AuthorList/Author")
  author_names <- vapply(authors, function(a) {
    fname <- xml_text(xml_find_first(a, "ForeName"))
    lname <- xml_text(xml_find_first(a, "LastName"))
    if (!nzchar(fname) || !nzchar(lname)) return(NA_character_)
    initials <- paste0(substr(unlist(strsplit(fname, " ")), 1, 1), collapse = "")
    paste0(lname, ", ", initials)
  }, character(1))
  author_str <- paste(na.omit(author_names), collapse = "; ")

  data.frame(
    PMID     = pmid,
    Title    = title,
    Abstract = abstract,
    Authors  = author_str,
    Journal  = journal,
    JournalAbbrev = journal_abbrev,
    Year     = year,
    stringsAsFactors = FALSE
  )
}

#' Fetch cardiovascular RCTs from PubMed.
#'
#' Uses `rentrez::entrez_search` + `entrez_fetch` to pull randomised controlled
#' trials from the target cardiology / general-medicine journals, filtered by a
#' cardiovascular MeSH / keyword scope.
#'
#' @param journals character vector of journal abbreviations. Defaults to the
#'   six target journals for the systematic CV-RCT pipeline.
#' @param start,end ISO-style `YYYY/MM/DD` date strings bounding `[PDAT]`.
#' @param out_csv path for the output CSV; `NULL` skips writing.
#' @param retmax hard upper bound on results (PubMed default is 20).
#' @param verbose logical, print progress.
#' @return data.frame of parsed article metadata.
#' @export
fetch_cv_rcts <- function(journals = c("N Engl J Med", "Lancet", "BMJ", "JAMA",
                                       "JAMA Cardiol", "Circulation"),
                          start    = "2020/01/01",
                          end      = "3000/12/31",
                          out_csv  = "data/pubmed_citations.csv",
                          retmax   = 5000,
                          verbose  = TRUE) {
  journal_clause <- .build_journal_clause(journals)
  cv_clause      <- .build_cv_clause()

  query <- paste(
    "(", journal_clause, ")",
    "AND", "(", cv_clause, ")",
    'AND "randomized controlled trial"[Publication Type]',
    'AND "humans"[MeSH Terms]',
    sprintf('AND (%s:%s[PDAT])', start, end),
    'NOT ("case reports"[Publication Type] OR "letter"[Publication Type])',
    'NOT ("systematic review"[Publication Type] OR "meta-analysis"[Publication Type])'
  )

  if (verbose) cat("PubMed query:\n", query, "\n\n", sep = "")

  init <- entrez_search(db = "pubmed", term = query, retmax = 0)
  total <- init$count
  if (verbose) cat("Total matching articles:", total, "\n")
  if (total == 0) {
    df <- data.frame(PMID = character(), Title = character(), Abstract = character(),
                     Authors = character(), Journal = character(),
                     JournalAbbrev = character(), Year = character(),
                     stringsAsFactors = FALSE)
    if (!is.null(out_csv)) {
      dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(df, out_csv, row.names = FALSE)
    }
    return(df)
  }

  n <- min(total, retmax)
  res <- entrez_search(db = "pubmed", term = query, retmax = n, sort = "pub date")
  pmids <- res$ids
  if (verbose) cat("Fetched", length(pmids), "PMIDs\n")

  # Fetch in batches of 200 to keep NCBI happy.
  batches <- split(pmids, ceiling(seq_along(pmids) / 200))
  parsed_all <- list()
  for (i in seq_along(batches)) {
    if (verbose) cat(sprintf("  batch %d/%d (%d records)\n",
                             i, length(batches), length(batches[[i]])))
    xml_raw <- entrez_fetch(db = "pubmed", id = batches[[i]],
                            rettype = "xml", parsed = FALSE)
    xml_doc <- read_xml(xml_raw)
    records <- xml_find_all(xml_doc, ".//PubmedArticle")
    parsed_all[[i]] <- do.call(rbind, lapply(records, .parse_article))
    Sys.sleep(0.34)  # be polite to NCBI (≤3 rps)
  }
  df <- do.call(rbind, parsed_all)

  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, out_csv, row.names = FALSE)
    if (verbose) cat("Saved metadata to", out_csv, "\n")
  }
  df
}
