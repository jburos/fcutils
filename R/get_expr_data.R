
#' List files containing featurecounts summaries by gene
#' @param featurecounts_dir path to directory containing featureCounts output files
#' @param pattern If provided, regex pattern for results files (e.g.: 'geneID$')
#' @param expect_length If provided, code will confirm the number of result files identified
#' @param sample_name_from_filename function to derive sample name from file name. Defaults to \code{basename}, but a useful alternative might be \code{dirname}.
#' @return named list of filepaths in directory matching regex, where name is extracted from the cleaned filename
#' @importFrom dplyr %>%
#' @import tidyverse
#' @export
get_featurecounts_files <- function(featurecounts_dir, pattern = NULL, expect_length = NULL, sample_name_from_filename = basename) {
  feature_counts <- dir(featurecounts_dir, full.names = T, include.dirs = T, recursive = T) %>%
    purrr::set_names(clean_file_names(sample_name_from_filename(.))) %>%
    unlist()
  if (!is.null(pattern))
    feature_counts <- 
      purrr::map(feature_counts, ~ purrr::keep(., stringr::str_detect, pattern = pattern))
  if (!is.null(expect_length))
    assertthat::assert_that(length(feature_counts) == expect_length,
                            msg = glue::glue('Incorrect number of featurecount files identified (expected {expect_length}).'))
  feature_counts
}

#' Clean file names by removing longest common substring from filename
#' @param x names to be cleaned
#' @return string of same dimension of x, with longest common substring removed
#' @export
clean_file_names <- function(x) {
  assertthat::assert_that(all(!duplicated(x)), msg = 'Error: file names are not unique. Check whether `sample_name_from_filename` function is correct.')
  common_part <- longest_common_substring(x)
  if (length(common_part) == 1)
    x <- stringr::str_remove(x, pattern = common_part)
  x
}

#' load data from featurecounts
#' borrowed from https://www.biostars.org/p/277316/#277350
#' @import purrr
#' @export
load_data_from_feature_counts <- function(filter = TRUE, files = NULL, featurecounts_dir, ...) {
  if (is.null(files))
    files <- get_featurecounts_files(featurecounts_dir = featurecounts_dir, ...)
  # list of tables (one for each file)
  l <- files %>% 
    purrr::map(~ read.table(., skip=1, header = T, stringsAsFactors = F)) %>%
    purrr::map(~ purrr::set_names(., c(names(.)[-1*ncol(.)], 'count')))
  if (!all(sapply(l, function(a) all(a$Geneid == l[[1]]$Geneid)))) 
    stop("Gene IDs (first column) differ between files.")
  # construct count matrix
  counts <- sapply(l, function(a) a$count)
  # colnames already pulled from names(files)
  rownames(counts) <- l[[1]][[1]] # rownames taken from first column
  # construct annotation data
  # could use `annotation <- getInBuiltAnnotation`, but prob better to get info from featurecounts output directly
  annotation <- l %>%
    purrr::map(tibble::rowid_to_column, var = 'rowid') %>%
    dplyr::bind_rows(., .id = 'sample_id') %>%
    dplyr::select(-sample_id, -count) %>%
    dplyr::distinct() %>%
    dplyr::arrange(rowid) %>%
    dplyr::select(-rowid)
  featurecounts <- tibble::lst(counts = counts, annotation = annotation)
  if (filter == TRUE) {
    featurecounts <- filter_feature_counts(featurecounts)
  }
  featurecounts
}

filter_feature_counts <- function(featurecounts) {
  keep <- filter_low_counts(featurecounts$counts)
  featurecounts$counts <- featurecounts$counts[keep,]
  featurecounts$annotation <- featurecounts$annotation[keep,]
  featurecounts
}

#' return gene annotation from ensembl ids
#' @export
#' @import biomaRt
get_gene_annotation <- function(ensembl_ids) {
  require("biomaRt")
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
  
  # remove suffixes from lookup ids
  ens <- ensembl_ids
  ensLookup <- gsub("\\.[0-9]*$", "", ens)
  
  # get annotation info
  annot <- biomaRt::getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)
  
  # add rows for ensembl ids that failed lookup
  missing_results <- ensLookup[!ensLookup %in% annot$ensembl_gene_id]
  dummy_annot <- tbl_df(list(ensembl_gene_id = missing_results))
  annot <- bind_rows(annot, dummy_annot)

  # add original_id back in as search item (in case modified by gsub, above)
  annot <- data.frame(
    ens[match(annot$ensembl_gene_id, ensLookup)],
    annot, stringsAsFactors = F)
  colnames(annot)[1] <- "original_id"
  
  # re-order annot to match original order of ensembl_ids
  annot <- annot[match(annot$ensembl_gene_id, ensembl_ids),]
  annot
}

#' get gene names from entrez or ensembl ids
#' @export
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db 
get_gene_names <- function(ids, column = 'SYMBOL', keytype=c("ENSEMBL", "ENTREZID"), ...) {
  keytype <- match.arg(keytype)
  if (keytype == 'ENSEMBL') {
    annot <- get_gene_annotation(ensembl_ids = ids)
    annot <- as.data.frame(annot)
    rownames(annot) <- annot$original_id
    return(annot[ids, 'external_gene_name'])
  }
  AnnotationDbi::mapIds(org.Hs.eg.db,
                        as.character(ids),
                        keytype=keytype,
                        column=column,
                        ...)
}

#' Summarise expression data as DGEList object
#' @export
get_expr_data <- function(using = c('DESeq2', 'edgeR'), sample_data = NULL, design = NULL) {
  # clean up arg inputs
  using <- match.arg(using)
  # get featurecount data
  fc <- load_data_from_feature_counts()
  # create result object
  if (using == 'DESeq2') {
    create_deseq2_result(fc = fc, sample_data = sample_data, design = design)
  } else {
    create_edger_result(fc = fc, sample_data = sample_data) # design ignored
  }
}

#' summarize expr data as DESeqDataSet
#' importFrom DESeq2 DESeqDataSetFromMatrix
#' @export
create_deseq2_result <- function(fc, sample_data, design = ~ cond) {
  # create DeSeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = fc$counts,
                                        colData = sample_data,
                                        design = design)
  dds
}

#' Summarise expr data as DGEList object
#' importFrom edgeR DGEList
#' @export
create_edger_result <- function(fc, sample_data) {
  dgel <- DGEList(fc$counts,
                  group = sample_data$cond,
                  samples = sample_data,
                  genes = fc$annotation)
  dgel
}
