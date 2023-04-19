#' Get GDC TCGA MAF files
#'
#' @param project a TCGA project
#'
#' @return maf file
#' @export

getsnvmaf <- function(project){

  if (!file.exists("output_snv")) {
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

  # 只要肿瘤样本
  #clin_snv <- clin[!clin$sample_type == "Solid Tissue Normal", ]

  # 只要snp文件中有的样本
  clin_snv <- clin[clin$submitter_id %in% index, ]

  # clin中没有Tumor_Sample_Barcode这一列，直接添加一列
  clin_snv$Tumor_Sample_Barcode <- clin_snv$submitter_id

  save(snv,clin_snv, file = paste0("output_snv/",project,"_maf_clin.rdata"))
}
