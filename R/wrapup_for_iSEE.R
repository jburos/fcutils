
#' Prepare summarized experiment from deseq2 result, for use in iSEE
#' borrowed from https://gist.github.com/federicomarini/4a543eebc7e7091d9169111f76d59de1
#' @export
#' @import DESeq2
#' @import SummarizedExperiment
wrapup_for_iSEE <- function(dds, res) {
  # dds to vst
  vst <- DESeq2::vst(dds)
  
  # initialize the container
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = List(
      counts = counts(dds),
      normcounts = counts(dds, normalized = TRUE),
      vst_counts = assay(vst)
    )
  )
  
  # adding colData, taken directly from the DESeqDataSet object
  colData(se) <- colData(dds)
  
  # extract contrast info
  this_contrast <- sub(".*p-value: (.*)", "\\1", mcols(res, use.names=TRUE)["pvalue","description"])
  
  # getting the rowData from the dds itself
  rdd <- rowData(dds)
  
  # modifying in advance the DESeqResults object
  res$log10_baseMean <- log10(res$baseMean)
  res$log10_pvalue <- -log10(res$pvalue)
  # and for the rowData
  rdd$log10_dispersion <- log10(rdd$dispersion)
  
  # adding rowData to se
  rowData(se)[[paste0("DESeq2_",gsub(" ","_",this_contrast))]] <- res
  
  # merging in the existing rowData slot
  rowData(se) <- cbind(rowData(se), rdd)
  
  return(se)
}
