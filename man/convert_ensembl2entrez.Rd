% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_expr_data.R
\name{convert_ensembl2entrez}
\alias{convert_ensembl2entrez}
\title{Given featurecounts matrix & annotation, roll
results up to the entrez gene id, combining
ensembl records mapping to the same gene}
\usage{
convert_ensembl2entrez(counts, annotation_df,
  original_id = "original_id", gene_id = "gene_id")
}
\arguments{
\item{counts}{counts matrix with rows per gene (named) & columns per sample (named)}

\item{annotation_df}{annotation df containing at least the original & gene ids}

\item{original_id}{character field name in the annotation df corresponding to row names in fc$counts matrix}

\item{gene_id}{character field name in the annotation df corresponding to the levels at which we want counts aggregated}
}
\value{
fc object with updated counts matrix & annotation df filtered to distinct gene_ids
}
\description{
Given featurecounts matrix & annotation, roll
results up to the entrez gene id, combining
ensembl records mapping to the same gene
}
