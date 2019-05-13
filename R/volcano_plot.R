
#' generate volcano plot from deseq2 result
#' @import DESeq2
#' @export
volcano_plot <- function(DESeq2_result,
                         alpha = 0.01, # Threshold on the p-value
                         min_fold_change = 2
                         ) {
  # par(mfrow=c(1,2))
  
  # Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
  DESeq2_result$sig <- -log10(DESeq2_result$padj)
  sum(is.infinite(DESeq2_result$sig))
  DESeq2_result[is.infinite(DESeq2_result$sig),"sig"] <- 350

  # Select genes with a defined p-value (DESeq2 assigns NA to some genes)
  genes.to.plot <- !is.na(DESeq2_result$pvalue)
  range(DESeq2_result[genes.to.plot, "log2FoldChange"])

  ## Volcano plot of adjusted p-values
  cols <- densCols(DESeq2_result$log2FoldChange, DESeq2_result$sig)
  cols[DESeq2_result$pvalue ==0] <- "purple"
  DESeq2_result$pch <- 19
  DESeq2_result$pch[DESeq2_result$pvalue ==0] <- 6
  plot(DESeq2_result$log2FoldChange, 
       DESeq2_result$sig, 
       col=cols, panel.first=grid(),
       main="Volcano plot", 
       xlab="Effect size: log2(fold-change)",
       ylab="-log10(adjusted p-value)",
       pch=DESeq2_result$pch, cex=0.4)
  abline(v=0)
  if (min_fold_change > 0)
    abline(v=c(-1*log2(min_fold_change),log2(min_fold_change)), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  ## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
  gn.selected <- abs(DESeq2_result$log2FoldChange) > min_fold_change & DESeq2_result$padj < alpha 
  text(DESeq2_result$log2FoldChange[gn.selected],
       -log10(DESeq2_result$padj)[gn.selected],
       lab=DESeq2_result$symbol[gn.selected ], cex=0.6)
}