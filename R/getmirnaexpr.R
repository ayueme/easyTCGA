#' Get GDC TCGA miRNA expression set
#'
#' @param project a TCGA project
#'
#' @return expression set
#' @export

getmirnaexpr <- function(project) {
  if (!file.exists("output_miRNA_expr")) {
    dir.create("output_miRNA_expr")
  }

  query <- TCGAbiolinks::GDCquery(project,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "miRNA Expression Quantification"
  )

  TCGAbiolinks::GDCdownload(query)

  TCGAbiolinks::GDCprepare(query, save = T,
                           save.filename = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))

  load(file = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))

  mirna_expr_counts <- data[, c(1, seq(2, ncol(data), 3))]
  colnames(mirna_expr_counts)[-1] <- substr(colnames(mirna_expr_counts)[-1], 12, 39)
  save(mirna_expr_counts, file = paste0("output_miRNA_expr/", project, "_mirna_expr_counts.rdata"))
  utils::write.csv(mirna_expr_counts,
            paste0("output_miRNA_expr/", project, "_mirna_expr_counts.csv"),
            quote = F,row.names = F)
  mirna_expr_rpm <- data[, c(1, seq(3, ncol(data), 3))]
  colnames(mirna_expr_rpm)[-1] <- substr(colnames(mirna_expr_rpm)[-1], 32, 59)
  save(mirna_expr_rpm, file = paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.rdata"))
  utils::write.csv(mirna_expr_rpm,
            paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.csv"),
            quote = F,row.names = F)
}
