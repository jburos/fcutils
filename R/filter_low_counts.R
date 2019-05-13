#' Filter out genes/transcripts with low expression
#' @param e expression data
#' @param cpm filter on cpm threshold default, otherwise use min_value
#' @param min_value if cpm is false, remove transcripts with rowsum <= this
#' @importFrom Biobase exprs
#' @importFrom edgeR DGEList
#' @export
filter_low_counts <- function(e, cpm=TRUE, min_value = 0) {
  if (cpm) {
    # this is the filtering recommendation from Gordon Smyth
    # https://support.bioconductor.org/p/85511/#85514
    L <- min(colSums(e)/1e6) # minimum library size (used to scale min theshold)
    dgel <- edgeR::DGEList(counts = e)
    # keep transcripts with expression in >=3 samples
    # since that is how many replicates we have
    # (using threshold to define "expression")
    return(rowSums(edgeR::cpm(dgel) > 10/L) >= 3)
  } else {
    return(rowSums(e) > min_value)
  }
}
