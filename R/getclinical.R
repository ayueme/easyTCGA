#' Get GDC TCGA clinical data
#'
#' @description
#' Get TCGA clinical data, including indexed clinical data, general information,
#' survival, tumor stage, drug, radiation, eta.
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
#' @return rdata files which are saved in the directory of "output_clinical".
#' - **TCGA-XXX_clinicalXML_admin.rdata**: admin data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_drug.rdata**: drug data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_followUp.rdata**: follow up data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_newTumorEvent.rdata**: new tumor event data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_patient.rdata**: general clinical data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_radiation.rdata**: radiation data extracted from XML files.
#' - **TCGA-XXX_clinicalXML_stageEvent.rdata**: pathologic stage data extracted from XML files.
#' - **TCGA-XXX_clinical_indexed.rdata**: indexed clinical data.
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
#' [getclinical()] can retrieve both the indexed and XML clinical data.
#'
#' @references Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
#' Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
#' Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
#' analysis of TCGA data.” Nucleic Acids Research.
#'
#' @export
#'

getclinical <- function(project){
  if (!dir.exists("output_clinical")){dir.create("output_clinical")}

  message("=> Querying data. \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "bcr xml"
  )

  message("\n=> Downloading data. \n")
  TCGAbiolinks::GDCdownload(query)

  message("\n=> Parsing clinical data.\n")
  # sapply(c("patient","drug","follow_up","radiation","admin","stage_event","new_tumor_event"),
  #        function(x){
  #          clinicalXML_patient <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = x)
  #          save(clinicalXML_patient, file = paste0("output_clinical/", project, "_XML_",x,".rdata"))
  #
  #        }
  # )

  clinicalXML_patient <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "patient")
  save(clinicalXML_patient, file = paste0("output_clinical/", project, "_clinicalXML_patient.rdata"))
  clinicalXML_drug <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "drug")
  save(clinicalXML_drug, file = paste0("output_clinical/", project, "_clinicalXML_drug.rdata"))
  clinicalXML_followUp <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "follow_up")
  save(clinicalXML_followUp, file = paste0("output_clinical/", project, "_clinicalXML_followUp.rdata"))
  clinicalXML_radiation <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "radiation")
  save(clinicalXML_radiation, file = paste0("output_clinical/", project, "_clinicalXML_radiation.rdata"))
  clinicalXML_admin <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "admin")
  save(clinicalXML_admin, file = paste0("output_clinical/", project, "_clinicalXML_admin.rdata"))
  clinicalXML_stageEvent <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "stage_event")
  save(clinicalXML_stageEvent, file = paste0("output_clinical/", project, "_clinicalXML_stageEvent.rdata"))
  clinicalXML_newTumorEvent <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = "new_tumor_event")
  save(clinicalXML_newTumorEvent, file = paste0("output_clinical/", project, "_clinicalXML_newTumorEvent.rdata"))

  # indexed clinical not available for TCGA-TGCT
  if(!project == "TCGA-TCGCT"){
    clinical_indexed <- TCGAbiolinks::GDCquery_clinic(project)
    save(clinical_indexed, file = paste0("output_clinical/", project, "_clinical_indexed.rdata"))
  }

  message("\n=> Successful.\n")
}




