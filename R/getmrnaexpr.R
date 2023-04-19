#' Get GDC TCGA mRNA and lncRNA expression set and clinical info
#'
#' @param project a TCGA project
#'
#' @return expression set and clinical info
#' @export

getmrnaexpr <- function(project) {
  if (!file.exists("output_mRNA_lncRNA_expr")) {
    dir.create("output_mRNA_lncRNA_expr")
  }
  print("Querying begins. Make sure your network has access to GDC TCGA!")

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  print("Downloading begins. Make sure your network has access to GDC TCGA!")

  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)

  print("Downloading ends. Preparing begins.")

  TCGAbiolinks::GDCprepare(query, save = T, save.filename = paste0("output_mRNA_lncRNA_expr/", project, "_expr.rdata"))

  load(file = paste0("output_mRNA_lncRNA_expr/", project, "_expr.rdata"))

  se <- data

  clin_info <- as.data.frame(SummarizedExperiment::colData(se))
  save(clin_info, file = paste0("output_mRNA_lncRNA_expr/", project, "_clinical.rdata"))
  utils::write.csv(clin_info,
            paste0("output_mRNA_lncRNA_expr/", project, "_clinical.csv"),
            quote = F,row.names = F)

  rowdata <- SummarizedExperiment::rowData(se)
  se_mrna <- se[rowdata$gene_type == "protein_coding", ]
  se_lnc <- se[rowdata$gene_type == "lncRNA", ]
  expr_counts_mrna <- SummarizedExperiment::assay(se_mrna, "unstranded")
  expr_tpm_mrna <- SummarizedExperiment::assay(se_mrna, "tpm_unstrand")
  expr_fpkm_mrna <- SummarizedExperiment::assay(se_mrna, "fpkm_unstrand")
  expr_counts_lnc <- SummarizedExperiment::assay(se_lnc, "unstranded")
  expr_tpm_lnc <- SummarizedExperiment::assay(se_lnc, "tpm_unstrand")
  expr_fpkm_lnc <- SummarizedExperiment::assay(se_lnc, "fpkm_unstrand")
  symbol_mrna <- SummarizedExperiment::rowData(se_mrna)$gene_name
  symbol_lnc <- SummarizedExperiment::rowData(se_lnc)$gene_name


  mrna_expr_counts <- cbind(data.frame(symbol_mrna), as.data.frame(expr_counts_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_counts[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_counts)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_counts <- tmp[,c(-1,-2)]
  save(mrna_expr_counts,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_counts.rdata"))
  utils::write.csv(mrna_expr_counts,
            paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_counts.csv"),
            quote = F)

  mrna_expr_tpm <- cbind(data.frame(symbol_mrna), as.data.frame(expr_tpm_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_tpm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_tpm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_tpm <- tmp[,c(-1,-2)]
  save(mrna_expr_tpm,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_tpm.rdata"))
  utils::write.csv(mrna_expr_tpm,
            paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_tpm.csv"),
            quote = F)

  mrna_expr_fpkm <- cbind(data.frame(symbol_mrna), as.data.frame(expr_fpkm_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_fpkm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_fpkm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_fpkm <- tmp[,c(-1,-2)]
  save(mrna_expr_fpkm,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_fpkm.rdata"))
  utils::write.csv(mrna_expr_fpkm,
            paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_fpkm.csv"),
            quote = F)


  lncrna_expr_counts <- cbind(data.frame(symbol_lnc), as.data.frame(expr_counts_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_counts[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_counts)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_counts <- tmp[,c(-1,-2)]
  save(lncrna_expr_counts,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_counts.rdata"))
  utils::write.csv(lncrna_expr_counts,
            paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_counts.csv"),
            quote = F)

  lncrna_expr_tpm <- cbind(data.frame(symbol_lnc), as.data.frame(expr_tpm_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_tpm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_tpm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_tpm <- tmp[,c(-1,-2)]
  save(lncrna_expr_tpm,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_tpm.rdata"))
  utils::write.csv(lncrna_expr_tpm,
            paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_tpm.csv"),
            quote = F)

  lncrna_expr_fpkm <- cbind(data.frame(symbol_lnc), as.data.frame(expr_fpkm_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_fpkm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_fpkm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_fpkm <- tmp[,c(-1,-2)]
  save(lncrna_expr_fpkm,
       file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_fpkm.rdata"))
  utils::write.csv(lncrna_expr_fpkm,
            paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_fpkm.csv"),
            quote = F)
}



