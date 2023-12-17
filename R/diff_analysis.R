#' Differential analysis
#'
#' @description This function can automatically do differential analysis by 4
#' methods: `DESeq2`, `egdeR`, `limma` and `wilcoxon test`.
#' @param exprset expression matrix prepared by [getmrnaexpr()]/[getmirnaexpr()]/
#' [getmrnaexpr_xena()], or your own expression matrix. If use your own
#' expression matrix, the argument "group" must be provided.
#' @param project characters used as part of file name.
#' @param is_count is the expression matrix type count? Default is TRUE. If TRUE,
#' the DEA will be done by `DESeq2`, `egdeR` and `limma`; if FALSE, the DEA will
#' be done by `limma` and `wilcoxon test`.
#' @param logFC_cut threshold of logFC used for filtering the DEGs, default is 0.
#' @param pvalue_cut p value cutpoint of the DEGs, default is 1.
#' @param adjpvalue_cut adjusted p value cutpoint of the DEGs, default is 1.
#' @param group a factor with two levels, specifying the group of samples, which
#' should be the same length as your sample size. If use your own expression
#' matrix, it must be provided.
#' @param save save results to local, default is FALSE.
#'
#' @return a list of differential analysis results.
#'
#' @references Love MI, Huber W, Anders S (2014). “Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.
#'
#' Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for
#' differential expression analysis of digital gene expression data.”
#' Bioinformatics, 26(1), 139-140.
#'
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies.” Nucleic Acids Research, 43(7), e47.
#'
#' @export

diff_analysis <- function(exprset,
                          project,
                          is_count = TRUE,
                          logFC_cut = 0,
                          pvalue_cut = 1,
                          adjpvalue_cut = 1,
                          group = NULL,
                          save = FALSE) {
  # check group
  if (!is.null(group)) {
    if (!is.factor(group)) stop("The type of group should be factor.")
    metadata <- data.frame(sample_id = colnames(exprset),group = group)
  } else {
    group <- ifelse(as.numeric(substr(colnames(exprset), 14, 15)) < 10, "tumor", "normal")
    group <- factor(group, levels = c("normal", "tumor"))
    metadata <- data.frame(sample_id = colnames(exprset),group = group)
  }
  # prepare DE analysis
  res_diff <- list()
  # check type
  if (is_count) {
    message("=> Running DESeq2")
    ## deseq2
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = exprset,colData = metadata,design = ~group)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, tidy = T)
    names(res)[1] <- "genesymbol"
    deg_deseq2 <- stats::na.omit(res)
    deg_deseq2 <- subset(deg_deseq2, abs(log2FoldChange)>logFC_cut & padj<adjpvalue_cut & pvalue<pvalue_cut)
    res_diff[[1]] <- deg_deseq2
    names(res_diff)[[1]] <- "deg_deseq2"
    ## limma voom
    message("=> Running limma voom")
    y <- edgeR::DGEList(counts = exprset, group = group)
    keep <- edgeR::filterByExpr(y, group = group)
    y <- y[keep, keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    design <- stats::model.matrix(~group)
    rownames(design) <- colnames(y)
    colnames(design) <- levels(y)
    v <- limma::voom(y, design, normalize = "quantile")
    fit <- limma::lmFit(v, design)
    fit2 <- limma::eBayes(fit)
    DEG2 <- limma::topTable(fit2, coef = 2, n = Inf)
    deg_limma <- stats::na.omit(DEG2)
    deg_limma$genesymbol <- rownames(deg_limma)
    deg_limma <- subset(deg_limma, abs(logFC)>logFC_cut & adj.P.Val<adjpvalue_cut & P.Value<pvalue_cut)
    res_diff[[2]] <- deg_limma
    names(res_diff)[[2]] <- "deg_limma"
    ## edger
    message("=> Running edgeR")
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    DEG <- edgeR::topTags(qlf, n = Inf)
    deg_edger <- as.data.frame(DEG)
    deg_edger$genesymbol <- rownames(deg_edger)
    deg_edger <- subset(deg_edger, abs(logFC)>logFC_cut & FDR<adjpvalue_cut & PValue<pvalue_cut)
    res_diff[[3]] <- deg_edger
    names(res_diff)[[3]] <- "deg_edger"
  } else {
    exprset <- limma::normalizeBetweenArrays(exprset)

    ex <- exprset
    qx <- as.numeric(stats::quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

    if (LogC) {
      # ex[which(ex <= 0)] <- NaN
      exprset <- log2(ex + 0.1)
      message("=> log2(x+0.1) transform finished")
    }else{message("=> log2 transform not needed")}

    ## limma
    message("=> Running limma")
    design <- stats::model.matrix(~group)
    fit <- limma::lmFit(exprset, design)
    fit <- limma::eBayes(fit)
    deg_limma <- limma::topTable(fit, coef = 2, number = Inf)
    deg_limma <- stats::na.omit(deg_limma)
    deg_limma$genesymbol <- rownames(deg_limma)
    deg_limma <- subset(deg_limma, abs(logFC)>logFC_cut & adj.P.Val<adjpvalue_cut & P.Value<pvalue_cut)
    res_diff[[1]] <- deg_limma
    names(res_diff)[[1]] <- "deg_limma"

    ## wilcoxon
    message("=> Running wilcoxon test")
    exprset <- t(exprset)
    wilcox_res <- apply(exprset, 2, function(x) {
      stats::wilcox.test(x ~ group, exact = F, correct = F)$p.value
    })

    deg_wilcoxon <- as.data.frame(wilcox_res)
    names(deg_wilcoxon) <- "pvalue"
    deg_wilcoxon$genesymbol <- rownames(deg_wilcoxon)
    res_diff[[2]] <- deg_wilcoxon
    names(res_diff)[[2]] <- "deg_wilcoxon"
  }

  if (save) {
    if (!dir.exists("output_diff")) {dir.create("output_diff")}
    save(res_diff, file = paste0("output_diff/", project, "_diff_results.rdata"))
  }
  message("=> Analysis done.")
  return(res_diff)
}

