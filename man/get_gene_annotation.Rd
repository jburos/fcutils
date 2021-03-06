% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_expr_data.R
\name{get_gene_annotation}
\alias{get_gene_annotation}
\title{return gene annotation from ensembl ids}
\usage{
get_gene_annotation(ensembl_ids, attributes = c("ensembl_gene_id",
  "gene_biotype", "external_gene_name", "entrezgene_id", "hgnc_symbol"),
  extra_attributes = c(), host = fcutils_options()$ensembl_host)
}
\description{
return gene annotation from ensembl ids
}
