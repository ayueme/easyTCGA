#' Differential analysis
#'
#' @description This function can automatically do differential analysis by 3
#'     ways: DESeq2, egdeR and limma.
#' @param exprset counts expression matrix prepared by getmrnaexpr() or
#'     getmirnaexpr(), or your own counts expression matrix. If use your own
#'     expression matrix, the argument "group" must be provided
#' @param project a character used as the part of file name
#' @param group a factor with two levels, specifying the group of samples,
#'     which should be the same length as your sample size. If use your own
#'     expression matrix, it must be provided
#' @param save save results to local, default is "FALSE"
#'
#' @return a list of differential analysis results from DESeq2, limma, edgeR
#' @export

diff_analysis <- function(exprset,
                          project,
                          group = NULL,
                          save = FALSE
){
  if(!is.null(group)){
    metadata <- data.frame(sample_id = colnames(exprset),
                           group = group
    )
  } else {
    group <- ifelse(substr(colnames(exprset),14,15)<10,"tumor","normal")
    group <- factor(group, levels = c("normal","tumor"))
    metadata <- data.frame(sample_id = colnames(exprset),
                           group = group
    )
  }
  res_diff <- list()
  ## deseq2
  cli::cli_alert_info("Running DESeq2")
  dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = exprset,
                                         colData = metadata,
                                         design = ~ group
  )
  dds <- DESeq2::DESeq(dds1)
  res <- DESeq2::results(dds, tidy = T)
  names(res)[1] <- "genesymbol"
  deg_deseq2 <- stats::na.omit(res)
  res_diff[[1]] <- deg_deseq2
  names(res_diff)[[1]] <- "deg_deseq2"
  ## limma voom
  cli::cli_alert_info("Running limma voom")
  y <- edgeR::DGEList(counts = exprset, group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  y <- y[keep,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~ group)
  rownames(design)<-colnames(y)
  colnames(design)<-levels(y)
  v <- limma::voom(y, design, normalize = "quantile")
  fit <- limma::lmFit(v, design)
  fit2 <- limma::eBayes(fit)
  DEG2 <- limma::topTable(fit2, coef=2, n=Inf)
  deg_limma <- stats::na.omit(DEG2)
  deg_limma$genesymbol <- rownames(deg_limma)
  res_diff[[2]] <- deg_limma
  names(res_diff)[[2]] <- "deg_limma"
  ## edger
  cli::cli_alert_info("Running edgeR")
  y <- edgeR::estimateDisp(y,design)
  fit <- edgeR::glmQLFit(y,design)
  qlf <- edgeR::glmQLFTest(fit,coef=2)
  DEG <- edgeR::topTags(qlf, n = Inf)
  deg_edger <- as.data.frame(DEG)
  deg_edger$genesymbol <- rownames(deg_edger)
  res_diff[[3]] <- deg_edger
  names(res_diff)[[3]] <- "deg_edger"

  if(save){
    if (!dir.exists("output_diff")) {
      dir.create("output_diff")
    }
    save(res_diff, file = paste0("output_diff/",project,"_deg_results.rdata"))
  }
  cli::cli_alert_success("Analysis done.")
  return(res_diff)
}

