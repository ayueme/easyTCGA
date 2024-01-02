#' Get GDC TCGA MAF files
#'
#' @description This function can automatically query, download and prepare the
#' TCGA maf files(masked somatic mutation) with the corresponding clinical
#' information. The output can be directly used by [maftools::read.maf()] function.
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
#' TCGA-COAD and TCGA-READ. Attention: there is no "matched_maf_and_clin" file
#' in this case.
#'
#' @return maf file and clinical information. The data are saved in the directory
#' of "output_snv".
#' - **TCGA-XXX_maf.rdata**: (masked somatic mutation)maf files
#' - **TCGA-XXX_clinical_indexed.rdata**: indexed clinical data
#' - **TCGA-XXX_matched_maf_and_clin.rdata**: matched maf file and clinical data,
#' the number and order of samples are exactly the same, so they can be directly
#' used by [maftools::read.maf()].
#'
#' @details
#' In GDC database the clinical data can be retrieved from different sources:
#' - indexed clinical: a refined clinical data that is created using the XML files
#' - XML files: original source of the data
#' - BCR Biotab: tsv files parsed from XML files
#'
#' more information about clinical data, please see: [clinical section](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html)
#' of the `TCGAbiolinks` vignette.
#'
#' @references Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
#' Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
#' Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
#' analysis of TCGA data.” Nucleic Acids Research.
#'
#' @export

getsnvmaf <- function(project){
  if (!dir.exists("output_snv")){dir.create("output_snv")}
  message("=> Querying data. \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open")
  message("\n=> Downloading data. \n")
  TCGAbiolinks::GDCdownload(query)

  if(length(project)<2){
    message("\n=> Preparing data.")
    TCGAbiolinks::GDCprepare(
      query, save = T,
      save.filename = paste0("output_snv/",project,"_maf.rdata"))
    clinical_indexed <- TCGAbiolinks::GDCquery_clinic(project = project,
                                                      type = "clinical")
    save(clinical_indexed,
         file = paste0("output_snv/",project,"_clinical_indexed.rdata"))
    load(file = paste0("output_snv/",project,"_maf.rdata"))
    #snv <- data
    data$Tumor_Sample_Barcode <- substr(data$Tumor_Sample_Barcode,1,12)
    index <- unique(data$Tumor_Sample_Barcode)
    clin_snv <- clinical_indexed[clinical_indexed$submitter_id %in% index, ]
    clin_snv$Tumor_Sample_Barcode <- clin_snv$submitter_id
    save(data,clin_snv,
         file = paste0("output_snv/",project,"_matched_maf_and_clin.rdata"))
    message("\n=> Successful.")
  }
  else{
    project <- paste(project,collapse = "_")
    message("\n=> Preparing data.")
    TCGAbiolinks::GDCprepare(
      query, save = T,
      save.filename = paste0("output_snv/",project,"_maf.rdata"))
    message("\n=> Successful.")
  }
}
