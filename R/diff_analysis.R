#' Differential analysis
#'
#' @param exprset expression set prepared by getmrnaexpr()
#' @param project a TCGA project
#' @param save save results to local, defautl is "FALSE"
#'
#' @return the differential analysis results
#' @export

diff_analysis <- function(exprset, project, save = FALSE){

  if (!file.exists("output_diff")) {
    dir.create("output_diff")
  }

  res_diffe <- list()
  ## deseq2
  group <- ifelse(substr(colnames(exprset),14,15)<10,"tumor","normal")
  group <- factor(group, levels = c("normal","tumor"))

  metadata <- data.frame(sample_id = colnames(exprset),
                         group = group
  )

  dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = exprset,
                                         colData = metadata,
                                         design = ~ group
  )
  dds <- DESeq2::DESeq(dds1)
  res <- DESeq2::results(dds, tidy = T)

  names(res)[1] <- "genesymbol"
  deg_deseq2 <- stats::na.omit(res)

  res_diffe[[1]] <- deg_deseq2
  names(res_diffe)[[1]] <- "deg_deseq2"

  ## limma voom

  y <- edgeR::DGEList(counts = exprset, group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  y <- y[keep,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ group)
  rownames(design)<-colnames(y)
  colnames(design)<-levels(y)

  v <- limma::voom(y, design, normalize = "quantile", plot = T)
  fit <- limma::lmFit(v, design)
  fit2 <- limma::eBayes(fit)

  DEG2 <- limma::topTable(fit2, coef=2, n=Inf)
  deg_limma <- stats::na.omit(DEG2)
  deg_limma$genesymbol <- rownames(deg_limma)

  res_diffe[[2]] <- deg_limma
  names(res_diffe)[[2]] <- "deg_limma"

  ## edger
  y <- edgeR::estimateDisp(y,design)

  fit <- edgeR::glmQLFit(y,design)
  qlf <- edgeR::glmQLFTest(fit,coef=2)

  DEG <- edgeR::topTags(qlf, n = Inf)
  deg_edger <- as.data.frame(DEG)
  deg_edger$genesymbol <- rownames(deg_edger)

  res_diffe[[3]] <- deg_edger
  names(res_diffe)[[3]] <- "deg_edger"

  return(res_diffe)

  if(save){
    save(res_diffe, file = paste0("output_diff/",project,"_deg_results.rdata"))
  }

}













