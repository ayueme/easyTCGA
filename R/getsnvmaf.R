#' Get GDC TCGA MAF files
#'
#' @description This function can automatically download and prepare the TCGA
#'     maf files(masked somatic mutation) with the corresponding clinical
#'     information. The output can be directly used by maftools::read.maf()
#'     function.
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
#' @return maf file and clinical information, which can be directly used by
#'     maftools::read.maf(). The data are saved in the directory of "output_snv"
#' @export

getsnvmaf <- function(project){
  if (!dir.exists("output_snv")) {
    dir.create("output_snv")
  }
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
  )
  TCGAbiolinks::GDCdownload(query)
  TCGAbiolinks::GDCprepare(query, save = T,
             save.filename = paste0("output_snv/",project,"_maf.rdata"))
  clin <- TCGAbiolinks::GDCquery_clinic(project = project, type = "clinical")
  save(clin,file = paste0("output_snv/",project,"_clin.rdata"))
  load(file = paste0("output_snv/",project,"_maf.rdata"))
  snv <- data
  snv$Tumor_Sample_Barcode <- substr(snv$Tumor_Sample_Barcode,1,12)
  index <- unique(snv$Tumor_Sample_Barcode)
  clin_snv <- clin[clin$submitter_id %in% index, ]
  clin_snv$Tumor_Sample_Barcode <- clin_snv$submitter_id
  save(snv,clin_snv, file = paste0("output_snv/",project,"_maf_clin.rdata"))
}
