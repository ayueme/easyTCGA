#' Get GDC TCGA miRNA expression matrix and clinical information
#'
#' @description This function provides a powerful workflow to query, download
#' and prepare the newest TCGA miRNA expression quantification data. All you
#' have to supply is a valid TCGA project name. It will automatically save two
#' types of expression matrix(count and rpm), both in rdata and csv formats, and
#' the corresponding clinical information.
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
#' TCGA-COAD and TCGA-READ. Attention: there is no "matched_clinical" file in
#' this case.
#'
#' @return miRNA expression matrix and clinical information. The data are saved
#'     in the directory of "output_miRNA_expr".
#' - **TCGA-XXX_miRNA.rdata**: original miRNA expression matrix, count and rpm
#'   file are extracted from this file.
#' - **TCGA-XXX_mirna_expr_count.rdata**: miRNA count expression matrix
#' - **TCGA-XXX_mirna_expr_rpm.rdata**: miRNA rpm expression matrix
#' - **TCGA-XXX_clinical_indexed.rdata**: indexed clinical data
#' - **TCGA-XXX_matched_clinical.rdata**: refined clinical data from indexed
#'   clinical data, so the number and order of samples are exactly the same as
#'   them in count and rpm expression matrix.
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

getmirnaexpr <- function(project) {
  if (!dir.exists("output_miRNA_expr")) {dir.create("output_miRNA_expr")}
  message("=> Querying data. \n")
  query <- TCGAbiolinks::GDCquery(project,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "miRNA Expression Quantification")
  message("\n=> Downloading data. \n")
  TCGAbiolinks::GDCdownload(query)
  if(length(project)>1){
    project <- paste(project,collapse = "_")
    message("\n=> Preparing miRNA data.")
    TCGAbiolinks::GDCprepare(
      query, save = T,
      save.filename = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))
    message("\n=> Preparing count and rpm.")
    load(file = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))
    #tmp <- data
    mirna_expr_count <- data[, c(1, seq(2, ncol(data), 3))]
    colnames(mirna_expr_count)[-1] <- substr(colnames(mirna_expr_count)[-1], 12, 39)
    rownames(mirna_expr_count) <- mirna_expr_count[,1]
    mirna_expr_count <- mirna_expr_count[,-1]
    save(mirna_expr_count,
         file = paste0("output_miRNA_expr/", project, "_mirna_expr_count.rdata"))
    utils::write.csv(mirna_expr_count,
                     paste0("output_miRNA_expr/", project, "_mirna_expr_count.csv"),
                     quote = F,row.names = T)
    mirna_expr_rpm <- data[, c(1, seq(3, ncol(data), 3))]
    colnames(mirna_expr_rpm)[-1] <- substr(colnames(mirna_expr_rpm)[-1], 32, 59)
    rownames(mirna_expr_rpm) <- mirna_expr_rpm[,1]
    mirna_expr_rpm <- mirna_expr_rpm[,-1]
    save(mirna_expr_rpm,
         file = paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.rdata"))
    utils::write.csv(mirna_expr_rpm,
                     paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.csv"),
                     quote = F,row.names = T)
    message("\n=> Successful.")
  }
  else{
    message("\n=> Preparing miRNA data.")
    TCGAbiolinks::GDCprepare(
      query, save = T,
      save.filename = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))
    message("\n=> Preparing count and rpm.")
    load(file = paste0("output_miRNA_expr/", project, "_miRNA.rdata"))
    data <- data
    mirna_expr_count <- data[, c(1, seq(2, ncol(data), 3))]
    colnames(mirna_expr_count)[-1] <- substr(colnames(mirna_expr_count)[-1], 12, 39)
    rownames(mirna_expr_count) <- mirna_expr_count[,1]
    mirna_expr_count <- mirna_expr_count[,-1]
    save(mirna_expr_count,
         file = paste0("output_miRNA_expr/", project, "_mirna_expr_count.rdata"))
    utils::write.csv(mirna_expr_count,
                     paste0("output_miRNA_expr/", project, "_mirna_expr_count.csv"),
                     quote = F,row.names = T)
    mirna_expr_rpm <- data[, c(1, seq(3, ncol(data), 3))]
    colnames(mirna_expr_rpm)[-1] <- substr(colnames(mirna_expr_rpm)[-1], 32, 59)
    rownames(mirna_expr_rpm) <- mirna_expr_rpm[,1]
    mirna_expr_rpm <- mirna_expr_rpm[,-1]
    save(mirna_expr_rpm,
         file = paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.rdata"))
    utils::write.csv(mirna_expr_rpm,
                     paste0("output_miRNA_expr/", project, "_mirna_expr_rpm.csv"),
                     quote = F,row.names = T)

    message("\n=> Preparing clinical data.")
    clinical_indexed <- TCGAbiolinks::GDCquery_clinic(project = project,
                                                      type = "clinical")
    save(clinical_indexed, # 这个要不要保存呢？
         file = paste0("output_miRNA_expr/",project,"_clinical_indexed.rdata"))
    clin_matched <- clinical_indexed[match(substr(colnames(mirna_expr_count),1,12),
                                           clinical_indexed$submitter_id),]
    save(clin_matched,
         file = paste0("output_miRNA_expr/",project,"_matched_clinical.rdata"))
    utils::write.csv(clin_matched,
                     paste0("output_miRNA_expr/",project,"_matched_clinical.csv"),
                     quote = F,row.names = F)
    message("\n=> Successful.")
  }

}

