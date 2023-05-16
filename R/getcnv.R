#' Get GDC TCGA copy number variation data
#'
#' @description This function can easily get the copy number variation data from
#'    GDC TCGA.
#'
#' @param project valid TCGA project name(s) from 33 TCGA projects
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
#' @return copy number variation data. The data are saved in the directory of
#'    "output_cnv".
#' @export
#'

getcnv <- function(project){
  if (!dir.exists("output_cnv")){dir.create("output_cnv")}
  cli::cli_alert_info("Querying begins. Make sure your network has access to GDC TCGA! \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",
    access = "open")
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)
  TCGAbiolinks::GDCprepare(query, save = T,save.filename = paste0("output_cnv/",project,"_CNV.rdata"))
}



#' Get GDC TCGA 450K DNA methylation beta value matirx
#'
#' @description This function can easily get the 450K DNA methylation beta value
#'    matrix and the corresponding clinical information from GDC TCGA.
#'
#' @param project valid TCGA project name(s) from 33 TCGA projects
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
#' @return 450k DNA methylation beta value matrix and clinical information and
#'    probe information. The data are saved in the directory of "output_methy".
#' @export
#'

getmethybeta <- function(project){
  if (!dir.exists("output_methy")){dir.create("output_methy")}
  cli::cli_alert_info("Querying begins. Make sure your network has access to GDC TCGA! \n")

  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450"
  )
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)
  TCGAbiolinks::GDCprepare(query,save = T,save.filename=paste0("output_methy/",project,"_methy_se.rdata"))

  load(file = paste0("output_methy/",project,"_methy_beta_se.rdata"))
  clin_info <- as.data.frame(SummarizedExperiment::colData(data))
  save(clin_info, file = paste0("output_methy/",project,"_clinical.rdata"))

  probe_info <- as.data.frame(SummarizedExperiment::rowData(data))
  save(probe_info, file = paste0("output_methy/",project,"_probe_info.rdata"))

  beta_expr <- SummarizedExperiment::assay(data)
  #beta_expr <- subset(beta_expr, subset = (rowSums(is.na(beta_expr)) == 0))
  pd <- data.frame(Sample_Name = colnames(beta_expr),
                   sample_type = ifelse(as.numeric(substr(colnames(beta_expr),14,15))<10,"tumor","normal"))
  save(beta_expr,pd, file = paste0("output_methy/",project,"_beta_expr_and_pd.rdata"))

}

