
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
  feature_counts <- dir(featurecounts_dir, full.names = T, include.dirs = T, recursive = T) %>% unlist()
  if (!is.null(pattern))
    feature_counts <- 
      purrr::map(feature_counts, ~ purrr::keep(., stringr::str_detect, pattern = pattern)) %>%
      purrr::compact() %>% 
      unlist()
  if (!is.null(sample_name_from_filename))
    feature_counts <- feature_counts %>% 
      purrr::set_names(clean_file_names(sample_name_from_filename(.)))
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

#' return gene annotation from ensembl (or other gene) ids
#' @export
#' @import biomaRt
get_gene_annotation <- function(lookup_ids, 
                                filter = "ensembl_gene_id",
                                attributes = c('gene_biotype', 'external_gene_name', 'entrezgene_id', 'hgnc_symbol', 'entrezgene_description', 'entrezgene_accession'),
                                extra_attributes = c(),
                                host = fcutils_options()$ensembl_host) {
  require("biomaRt")
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart)
  
  if (length(extra_attributes) > 0)
    attributes <- unique(c(attributes, extra_attributes))
  attributes <- unique(c(attributes, filter))

  # remove suffixes from lookup ids
  if (filter == 'ensembl_gene_id') {
    # extract numeric identifier from EN* string
    lookup_ids <- gsub("\\.[0-9]*$", "", lookup_ids)
  }
  
  # get annotation info
  annot <- biomaRt::getBM(
    mart=mart,
    attributes=attributes,
    filter=filter,
    values=lookup_ids,
    uniqueRows=TRUE)
  
  # add rows for ensembl ids that failed lookup
  missing_results <- lookup_ids[!lookup_ids %in% annot[[filter]]]
  if (length(missing_results) > 0) {
    futile.logger::flog.warn(glue::glue('{length(missing_results)} ids not found in biomaRt. Dummy records (with all NA values) will be added to the annotation df.'))
    dummy_annot <- tbl_df(list(id = missing_results))
    names(dummy_annot) <- filter
    annot <- bind_rows(annot, dummy_annot)
  }
  
  # add original_id back in as search item (in case modified by gsub, above)
  annot <- annot %>%
    dplyr::mutate(original_id = !!rlang::sym(filter)) %>%
    as.data.frame()
  # check for duplicates
  if ((n_dups <- annot %>% dplyr::add_count(original_id) %>% filter(n>1) %>% distinct(original_id) %>% nrow()) > 1) {
    futile.logger::flog.warn(glue::glue('{n_dups} records with duplicates by {filter} in annotation results. These will be kept.'))
  }

  annot
}

#' get gene names from entrez or ensembl ids
#' @export
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db 
get_gene_names <- function(ids, column = 'SYMBOL', keytype=c("ENSEMBL", "ENTREZID"), host = fcutils_options()$ensembl_host, ...) {
  keytype <- match.arg(keytype)
  if (keytype == 'ENSEMBL') {
    annot <- get_gene_annotation(lookup_ids = ids, filter = 'ensembl_gene_id', host = host)
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

#' Given featurecounts matrix & annotation, roll
#' results up to the entrez gene id, combining 
#' ensembl records mapping to the same gene
#' @param counts counts matrix with rows per gene (named) & columns per sample (named)
#' @param annotation_df annotation df containing at least the original & gene ids
#' @param original_id character field name in the annotation df corresponding to row names in fc$counts matrix
#' @param gene_id character field name in the annotation df corresponding to the levels at which we want counts aggregated
#' @return fc object with updated counts matrix & annotation df filtered to distinct gene_ids
#' @export
convert_ensembl2entrez <- function(counts, annotation_df, original_id = 'original_id', gene_id = 'gene_id') {
  if (missing(annotation_df)) {
    annotation_df <- get_gene_annotation(rownames(counts))
    original_id <- 'original_id'
    gene_id <- 'hgnc_symbol'
    if (length(not_found <- rownames(counts)[!rownames(counts) %in% annotation_df[[original_id]]])>0) {
      futile.logger::flog.warn(glue::glue('{length(not_found)} gene ids could not be located in biomaRt: {glue::glue_collapse(head(not_found, n = 10), sep = ", ")}'))
    }
  }
  if (any(is.na(annotation_df[[original_id]])) || any(is.na(annotation_df[[gene_id]]))) {
    annotation_df_new <- annotation_df %>%
      dplyr::select(one_of(original_id, gene_id)) %>%
      na.omit()
    futile.logger::flog.info(glue::glue('Dropping records with missing {original_id} or {gene_id}; {nrow(annotation_df_new)}/{nrow(annotation_df)} remaining.'))
    annotation_df <- annotation_df_new
  }
  new_genes <- unique(annotation_df[[gene_id]])
  new_counts <- matrix(NA_integer_,
                       nrow = length(new_genes), ncol = ncol(counts),
                       dimnames = list(gene = new_genes,
                                       sample = colnames(counts)))
  for (gene in new_genes) {
    ensembl_ids <- annotation_df[annotation_df[[gene_id]] == gene, original_id]
    if (length(ensembl_ids) == 1) {
      new_counts[gene,] <- counts[ensembl_ids,]
    } else {
      subset <- counts[ensembl_ids,]
      new_counts[gene,] <- colSums(subset, na.rm = T)
    }
  }
  new_counts
}