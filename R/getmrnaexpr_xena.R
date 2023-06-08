#' Get mRNA and lncRNA expression matrix from XENA file
#'
#' @param expr_file gene expression RNAseq file downloaded from XENA(gdchub)
#' @param clin_file phenotype file downloaded from XENA(gdchub)
#' @return mRNA, lncRNA expression matrix and clinical info. The data are saved
#'    under the directory of "output_mRNA_expr_xena".
#' @export

getmrnaexpr_xena <- function(expr_file, clin_file = NULL){

  message("=> Make sure your files are downloaded from https://gdc.xenahubs.net.")
  if(!dir.exists("output_mrna_expr_xena")){dir.create("output_mrna_expr_xena")}

  exprset <- utils::read.table(expr_file,sep = "\t",header = T,check.names = F)
  # gtf v22, keep the same as XENA
  # mRNA
  gencodev22_mrna <- subset(genecodev22,
                            type=="gene"&gene_type=="protein_coding",
                            select = c("gene_id","gene_name")
  )
  expr_mrna <- merge(gencodev22_mrna,exprset, by.x = "gene_id",by.y = "Ensembl_ID")
  expr_mrna <- expr_mrna[,-1]
  rowm <- as.data.frame(rowMeans(expr_mrna[, -1]))
  names(rowm) <- "rmea"
  expr_mrna <- cbind(rowm, expr_mrna)
  expr_mrna <- expr_mrna[order(expr_mrna$rmea,decreasing=T),]
  expr_mrna <- expr_mrna[!duplicated(expr_mrna$gene_name),]
  rownames(expr_mrna) <- expr_mrna[,2]
  expr_mrna <- expr_mrna[,c(-1,-2)]

  # lncRNA
  # http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
  gencodev22_lncrna <- subset(genecodev22,
                              gene_type %in% c("3prime_overlapping_ncrna",
                                               "non_coding","antisense","lincRNA",
                                               "sense_intronic","sense_overlapping",
                                               "macro_lncRNA"),
                              select = c("gene_id","gene_name")
  )
  expr_lncrna <- merge(gencodev22_lncrna,exprset, by.x = "gene_id",by.y = "Ensembl_ID")
  expr_lncrna <- expr_lncrna[,-1]
  table(duplicated(expr_lncrna$gene_name))
  rowm <- as.data.frame(rowMeans(expr_lncrna[, -1]))
  names(rowm) <- "rmea"
  expr_lncrna <- cbind(rowm, expr_lncrna)
  expr_lncrna <- expr_lncrna[order(expr_lncrna$rmea,decreasing=T),]
  expr_lncrna <- expr_lncrna[!duplicated(expr_lncrna$gene_name),]
  rownames(expr_lncrna) <- expr_lncrna[,2]
  expr_lncrna <- expr_lncrna[,c(-1,-2)]

  if(grepl("counts",expr_file)){

    # check type
    expr_mrna <- 2^expr_mrna-1
    mrna_expr_count_xena <- floor(expr_mrna)

    expr_lncrna <- 2^expr_lncrna-1
    lncrna_expr_count_xena <- floor(expr_lncrna)

    save(mrna_expr_count_xena, file = "output_mrna_expr_xena/mrna_expr_count_xena.rdata")
    save(lncrna_expr_count_xena, file = "output_mrna_expr_xena/lncrna_expr_count_xena.rdata")
  }else{
    mrna_expr_fpkm_xena <- expr_lncrna
    lncrna_expr_fpkm_xena <- expr_lncrna
    save(mrna_expr_fpkm_xena, file = "output_mrna_expr_xena/mrna_expr_fpkm_xena.rdata")
    save(lncrna_expr_fpkm_xena, file = "output_mrna_expr_xena/lncrna_expr_fpkm_xena.rdata")
  }

  ## prepare clinical info from gdchub
  if(!is.null(clin_file)){

    # phenotype
    if(grepl("phenotype",clin_file)){
      clin_xena <- data.table::fread(clin_file,data.table = F)
      index <- intersect(colnames(expr_mrna),clin_xena$submitter_id.samples)
      expr_mrna <- expr_mrna[,index]
      clin_xena <- clin_xena[clin_xena$submitter_id.samples %in% index,]
      clin_xena <- clin_xena[match(colnames(expr_mrna),clin_xena$submitter_id.samples),]
      save(clin_xena, expr_mrna,file = "output_mrna_expr_xena/mrna_and_clin_xena.rdata")

    }else{
      # survival
      clin_xena <- data.table::fread(clin_file,data.table = F)
      index <- intersect(colnames(expr_mrna),clin_xena$sample)
      expr_mrna <- expr_mrna[,index]
      clin_xena <- clin_xena[clin_xena$sample %in% index,]
      clin_xena <- clin_xena[match(colnames(expr_mrna),clin_xena$sample),]
      save(clin_xena, expr_mrna,file = "output_mrna_expr_xena/mrna_and_clin_xena.rdata")
    }

  }

}
