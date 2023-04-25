#' Get GDC TCGA miRNA expression matrix
#'
#' @description This function provides a powerful workflow to query, download
#'     and prepare the newest TCGA miRNA expression quantification data. All you
#'     have to supply is a valid TCGA project name. It will automatically save
#'     two types of expression matrix(counts and rpm), both in rdata and csv
#'     formats.
#' @param project one of 33 TCGA projects
#' \itemize{
#' \item{ TCGA-ACC }
#' \item{ TCGA-BLCA }
#' \item{ TCGA-BRCA }
#' \item{ TCGA-CESC }
#' \item{ TCGA-CHOL }
#' \item{ TCGA-COAD }
#' \item{ TCGA-DLBC }
#' \item{ TCGA-ESCA }
#' \item{ TCGA-GBM }
#' \item{ TCGA-HNSC }
#' \item{ TCGA-KICH }
#' \item{ TCGA-KIRC }
#' \item{ TCGA-KIRP }
#' \item{ TCGA-LAML }
#' \item{ TCGA-LGG }
#' \item{ TCGA-LIHC }
#' \item{ TCGA-LUAD }
#' \item{ TCGA-LUSC }
#' \item{ TCGA-MESO }
#' \item{ TCGA-OV }
#' \item{ TCGA-PAAD }
#' \item{ TCGA-PCPG }
#' \item{ TCGA-PRAD }
#' \item{ TCGA-READ }
#' \item{ TCGA-SARC }
#' \item{ TCGA-SKCM }
#' \item{ TCGA-STAD }
#' \item{ TCGA-TGCT }
#' \item{ TCGA-THCA }
#' \item{ TCGA-THYM }
#' \item{ TCGA-UCEC }
#' \item{ TCGA-UCS }
#' \item{ TCGA-UVM }
#' }
#'
#' @return miRNA expression matrix. The data are saved in the directory of
#'     "output_miRNA_expr".
#' @export

getmirnaexpr <- function(project) {
  if (!dir.exists("output_miRNA_expr")) {
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
  rownames(mirna_expr_counts) <- mirna_expr_counts[,1]
  mirna_expr_counts <- mirna_expr_counts[,-1]
  save(mirna_expr_counts, file = paste0("output_miRNA_expr/", project, "_mirna_expr_counts.rdata"))
  utils::write.csv(mirna_expr_counts,
                   paste0("output_miRNA_expr/", project, "_mirna_expr_counts.csv"),
                   quote = F,row.names = T)
  mirna_expr_rpm <- data[, c(1, seq(3, ncol(data), 3))]
  colnames(mirna_expr_rpm)[-1] <- substr(colnames(mirna_expr_rpm)[-1], 32, 59)
  rownames(mirna_expr_rpm) <- mirna_expr_rpm[,1]
  mirna_expr_rpm <- mirna_expr_rpm[,-1]
  save(mirna_expr_rpm, file = paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.rdata"))
  utils::write.csv(mirna_expr_rpm,
                   paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.csv"),
                   quote = F,row.names = T)
}
