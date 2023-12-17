#' Get GDC TCGA mRNA/lncRNA expression matrix and clinical information
#'
#' @description This function provides a powerful workflow to query, download
#' and prepare the newest TCGA gene expression quantification data and the
#' clinical information. All you have to supply is a valid TCGA project name. It
#' can automatically save six types of expression matrix(mRNAcount/tpm/fpkm,
#' lncRNA count/tpm/fpkm) and the corresponding clinical information, both in
#' rdata and csv formats.
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
#' @return rdata and csv files which are saved in the directory of "output_mRNA_lncRNA_expr".
#' - **TCGA-XXX_SummarizedExperiment.rdata**: SummarizedExperiment object,
#'   all the other files are extracted from this object.
#' - **TCGA-XXX_clinicalSE.rdata**: indexed clinical information extracted
#'   from the SummarizedExperiment object.
#' - **TCGA-XXX_gene_info.rdata**: gene information, including HGNC Gene
#'   Symbol, Ensembl ID, gene type, eta.
#' - **TCGA-XXX_mrna_expr_count.rdata**: mRNA count expression matrix
#' - **TCGA-XXX_mrna_expr_tpm.rdata**: mRNA tpm expression matrix
#' - **TCGA-XXX_mrna_expr_fpkm.rdata**: mRNA fpkm expression matrix
#' - **TCGA-XXX_lncrna_expr_count.rdata**: lncRNA count expression matrix
#' - **TCGA-XXX_lncrna_expr_tpm.rdata**: lncRNA tpm expression matrix
#' - **TCGA-XXX_lncrna_expr_fpkm.rdata**: lncRNA fpkm expression matrix
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

getmrnaexpr <- function(project) {
  if (!dir.exists("output_mRNA_lncRNA_expr")){dir.create("output_mRNA_lncRNA_expr")}
  message("=> Querying begins. \n")
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - count"
  )
  message("\n=> Downloading begins. \n")
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 100)
  message("\n=> Preparing SummarizedExperiment.\n")
  if(length(project)>1){project <- paste(project,collapse = "_")}
  TCGAbiolinks::GDCprepare(query, save = T, save.filename = paste0("output_mRNA_lncRNA_expr/", project, "_SummarizedExperiment.rdata"))
  message("\n=> Preparing mRNA and lncRNA.")
  load(file = paste0("output_mRNA_lncRNA_expr/", project, "_SummarizedExperiment.rdata"))
  #se <- data
  clinicalSE <- as.data.frame(SummarizedExperiment::colData(data)) # 临床信息
  save(clinicalSE, file = paste0("output_mRNA_lncRNA_expr/", project, "_clinicalSE.rdata"))
  clinicalSE <- apply(clinicalSE,2,as.character)
  utils::write.csv(clinicalSE,paste0("output_mRNA_lncRNA_expr/", project, "_clinicalSE.csv"),row.names = F)

  rowdata <- SummarizedExperiment::rowData(data)
  # 基因名
  gene_info <- as.data.frame(rowdata@listData)
  save(gene_info,file = paste0("output_mRNA_lncRNA_expr/", project, "_gene_info.rdata"))
  utils::write.csv(gene_info,paste0("output_mRNA_lncRNA_expr/", project, "_gene_info.csv"),quote = F,row.names = F)

  se_mrna <- data[rowdata$gene_type == "protein_coding", ]
  se_lnc <- data[rowdata$gene_type == "lncRNA", ]
  expr_count_mrna <- SummarizedExperiment::assay(se_mrna, "unstranded")
  expr_tpm_mrna <- SummarizedExperiment::assay(se_mrna, "tpm_unstrand")
  expr_fpkm_mrna <- SummarizedExperiment::assay(se_mrna, "fpkm_unstrand")
  expr_count_lnc <- SummarizedExperiment::assay(se_lnc, "unstranded")
  expr_tpm_lnc <- SummarizedExperiment::assay(se_lnc, "tpm_unstrand")
  expr_fpkm_lnc <- SummarizedExperiment::assay(se_lnc, "fpkm_unstrand")
  symbol_mrna <- SummarizedExperiment::rowData(se_mrna)$gene_name
  symbol_lnc <- SummarizedExperiment::rowData(se_lnc)$gene_name

  mrna_expr_count <- cbind(data.frame(symbol_mrna), as.data.frame(expr_count_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_count[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_count)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_count <- tmp[,c(-1,-2)]
  save(mrna_expr_count,file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_count.rdata"))
  utils::write.csv(mrna_expr_count,paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_count.csv"),quote = F)

  mrna_expr_tpm <- cbind(data.frame(symbol_mrna), as.data.frame(expr_tpm_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_tpm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_tpm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_tpm <- tmp[,c(-1,-2)]
  save(mrna_expr_tpm,file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_tpm.rdata"))
  utils::write.csv(mrna_expr_tpm,paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_tpm.csv"),quote = F)

  mrna_expr_fpkm <- cbind(data.frame(symbol_mrna), as.data.frame(expr_fpkm_mrna))
  rowm <- as.data.frame(rowMeans(mrna_expr_fpkm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, mrna_expr_fpkm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_mrna),]
  rownames(tmp) <- tmp[,2]
  mrna_expr_fpkm <- tmp[,c(-1,-2)]
  save(mrna_expr_fpkm,file = paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_fpkm.rdata"))
  utils::write.csv(mrna_expr_fpkm,paste0("output_mRNA_lncRNA_expr/", project, "_mrna_expr_fpkm.csv"),quote = F)

  lncrna_expr_count <- cbind(data.frame(symbol_lnc), as.data.frame(expr_count_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_count[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_count)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_count <- tmp[,c(-1,-2)]
  save(lncrna_expr_count,file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_count.rdata"))
  utils::write.csv(lncrna_expr_count,paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_count.csv"),quote = F)

  lncrna_expr_tpm <- cbind(data.frame(symbol_lnc), as.data.frame(expr_tpm_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_tpm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_tpm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_tpm <- tmp[,c(-1,-2)]
  save(lncrna_expr_tpm,file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_tpm.rdata"))
  utils::write.csv(lncrna_expr_tpm,paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_tpm.csv"),quote = F)

  lncrna_expr_fpkm <- cbind(data.frame(symbol_lnc), as.data.frame(expr_fpkm_lnc))
  rowm <- as.data.frame(rowMeans(lncrna_expr_fpkm[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, lncrna_expr_fpkm)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$symbol_lnc),]
  rownames(tmp) <- tmp[,2]
  lncrna_expr_fpkm <- tmp[,c(-1,-2)]
  save(lncrna_expr_fpkm,file = paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_fpkm.rdata"))
  utils::write.csv(lncrna_expr_fpkm,paste0("output_mRNA_lncRNA_expr/", project, "_lncrna_expr_fpkm.csv"),quote = F)
  message("\n=> Successful.")
}



