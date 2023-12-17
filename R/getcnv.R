#' Get GDC TCGA copy number variation data
#'
#' @description This function can automatically query, download and prepare the
#' copy number variation data(Masked Copy Number Segment) from GDC TCGA.
#'
#' @param project valid TCGA project name from 33 TCGA projects:
#' - TCGA-ACC
#' - TCGA-BLCA
#' - TCGA-BRCA
#' - TCGA-CESC
#' - TCGA-CHOL
#' - TCGA-COAD
#' - TCGA-DLBC
#' - TCGA-ESCA
#' - TCGA-GBM
#' - TCGA-HNSC
#' - TCGA-KICH
#' - TCGA-KIRC
#' - TCGA-KIRP
#' - TCGA-LAML
#' - TCGA-LGG
#' - TCGA-LIHC
#' - TCGA-LUAD
#' - TCGA-LUSC
#' - TCGA-MESO
#' - TCGA-OV
#' - TCGA-PAAD
#' - TCGA-PCPG
#' - TCGA-PRAD
#' - TCGA-READ
#' - TCGA-SARC
#' - TCGA-SKCM
#' - TCGA-STAD
#' - TCGA-TGCT
#' - TCGA-THCA
#' - TCGA-THYM
#' - TCGA-UCEC
#' - TCGA-UCS
#' - TCGA-UVM
#'
#' @return copy number variation data. The data are saved in the directory of
#' "output_cnv".
#'
#' @references Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
#' Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
#' Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
#' analysis of TCGA data.” Nucleic Acids Research.
#'
#' @export
#'

getcnv <- function(project){
  if (!dir.exists("output_cnv")){dir.create("output_cnv")}
  message("=> Querying begins. \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",
    access = "open")
  message("=> Downloading begins. \n")
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)
  if(length(project)>1){project <- paste(project,collapse = "_")}
  message("\n=> Preparing begins.")
  TCGAbiolinks::GDCprepare(query, save = T,save.filename = paste0("output_cnv/",project,"_CNV.rdata"))
  message("\n=> Successful.")
}



#' Get GDC TCGA 450K DNA methylation beta value matirx
#'
#' @description This function can automatically query, download and prepare the
#' 450K DNA methylation beta value matrix and the corresponding clinical
#' information from GDC TCGA.
#'
#' @param project valid TCGA project name(s) from 33 TCGA projects:
#' - TCGA-ACC
#' - TCGA-BLCA
#' - TCGA-BRCA
#' - TCGA-CESC
#' - TCGA-CHOL
#' - TCGA-COAD
#' - TCGA-DLBC
#' - TCGA-ESCA
#' - TCGA-GBM
#' - TCGA-HNSC
#' - TCGA-KICH
#' - TCGA-KIRC
#' - TCGA-KIRP
#' - TCGA-LAML
#' - TCGA-LGG
#' - TCGA-LIHC
#' - TCGA-LUAD
#' - TCGA-LUSC
#' - TCGA-MESO
#' - TCGA-OV
#' - TCGA-PAAD
#' - TCGA-PCPG
#' - TCGA-PRAD
#' - TCGA-READ
#' - TCGA-SARC
#' - TCGA-SKCM
#' - TCGA-STAD
#' - TCGA-TGCT
#' - TCGA-THCA
#' - TCGA-THYM
#' - TCGA-UCEC
#' - TCGA-UCS
#' - TCGA-UVM
#'
#' If you provide more than one TCGA project names, it will combine the data.
#' This is useful when you try to combine different cancer types, such as
#' TCGA-COAD and TCGA-READ.
#'
#' @return rdata files which are saved in the directory of "output_methy".
#' - **TCGA-XXX_methy_beta_SummarizedExperiment.rdata**: SummarizedExperiment
#'   object, all the other files are extracted from this object.
#' - **TCGA-XXX_clinicalSE.rdata**: indexed clinical information extracted
#'   from the SummarizedExperiment object.
#' - **TCGA-XXX_probe_info.rdata**: probe information
#' - **TCGA-XXX_methy_beta_expr.rdata**: DNA methylation beta value matrix
#' - **TCGA-XXX_pd.rdata**: a simple pd file including sample name and sample
#'   type(normal or tumor), which can be used in `ChAMP`.
#'
#' @details
#' In GDC database the clinical data can be retrieved from different sources:
#' - indexed clinical: a refined clinical data that is created using the XML files
#' - XML files: original source of the data
#' - BCR Biotab: tsv files parsed from XML files
#'
#' the clinical data extracted from the SummarizedExperiment object is indexed
#' clinical, and subtype information from marker TCGA papers are added, these
#' information have prefix "paper_". more information about clinical data,
#' please see: [clinical section](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html) of the `TCGAbiolinks` vignette.
#'
#' @references Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
#' Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
#' Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
#' analysis of TCGA data.” Nucleic Acids Research.
#'
#' @export
#'

getmethybeta <- function(project){
  if (!dir.exists("output_methy")){dir.create("output_methy")}
  message("=> Querying begins. \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450"
  )
  message("=> Downloading begins. \n")
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)
  if(length(project)>1){project <- paste(project,collapse = "_")}
  message("\n=> Preparing begins.")
  TCGAbiolinks::GDCprepare(query,save = T,save.filename=paste0("output_methy/",project,"_methy_beta_SummarizedExperiment.rdata"))
  load(file = paste0("output_methy/",project,"_methy_beta_SummarizedExperiment.rdata"))
  clinicalSE <- as.data.frame(SummarizedExperiment::colData(data))
  save(clinicalSE, file = paste0("output_methy/",project,"_clinicalSE.rdata"))
  probe_info <- as.data.frame(SummarizedExperiment::rowData(data))
  save(probe_info, file = paste0("output_methy/",project,"_probe_info.rdata"))
  beta_expr <- SummarizedExperiment::assay(data)
  save(beta_expr, file = paste0("output_methy/",project,"_methy_beta_expr.rdata"))
  pd <- data.frame(Sample_Name = colnames(beta_expr),
                   sample_type = ifelse(as.numeric(substr(colnames(beta_expr),14,15))<10,"tumor","normal"))
  save(pd, file = paste0("output_methy/",project,"_pd.rdata"))
  message("\n=> Successful.")
}


