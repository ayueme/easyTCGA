#' Get mRNA and lncRNA expression matrix from XENA file
#'
#' @param expr_file gene expression RNAseq file downloaded from XENA(gdchub)
#' @param clin_file phenotype file downloaded from XENA(gdchub)
#' @return mRNA, lncRNA expression matrix and clinical info. The data are saved
#'    under the directory of "output_mRNA_expr_xena".
#' @export

getmrnaexpr_xena <- function(expr_file, clin_file = NULL){

  if(!dir.exists("output_mrna_expr_xena")){dir.create("output_mrna_expr_xena")}

  exprset <- utils::read.table(expr_file,sep = "\t",header = T,check.names = F)
  # gtf v22, keep the same as XENA
  # mRNA
  gencodev22_mrna <- subset(genecodev22,
                            type=="gene"&gene_type=="protein_coding",
                            select = c("gene_id","gene_name")
  )
  tmp <- merge(gencodev22_mrna,exprset, by.x = "gene_id",by.y = "Ensembl_ID")
  tmp <- tmp[,-1]
  rowm <- as.data.frame(rowMeans(tmp[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, tmp)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$gene_name),]
  rownames(tmp) <- tmp[,2]
  tmp <- tmp[,c(-1,-2)]

  if(grepl("counts",expr_file)){

    tmp <- 2^tmp-1
    mrna_expr_xena <- floor(tmp)
  }else{mrna_expr_xena <- tmp}

  save(mrna_expr_xena, file = "output_mrna_expr_xena/mrna_expr_xena.rdata")
  # lncRNA
  # http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
  gencodev22_lncrna <- subset(genecodev22,
                              gene_type %in% c("3prime_overlapping_ncrna",
                                               "non_coding","antisense","lincRNA",
                                               "sense_intronic","sense_overlapping",
                                               "macro_lncRNA"),
                              select = c("gene_id","gene_name")
  )
  tmp <- merge(gencodev22_lncrna,exprset, by.x = "gene_id",by.y = "Ensembl_ID")
  tmp <- tmp[,-1]
  table(duplicated(tmp$gene_name))
  rowm <- as.data.frame(rowMeans(tmp[, -1]))
  names(rowm) <- "rmea"
  tmp <- cbind(rowm, tmp)
  tmp <- tmp[order(tmp$rmea,decreasing=T),]
  tmp <- tmp[!duplicated(tmp$gene_name),]
  rownames(tmp) <- tmp[,2]
  tmp <- tmp[,c(-1,-2)]

  if(grepl("counts",expr_file)){

    tmp <- 2^tmp-1
    lncrna_expr_xena <- floor(tmp)
  }else{lncrna_expr_xena <- tmp}


  save(lncrna_expr_xena, file = "output_mrna_expr_xena/lncrna_expr_xena.rdata")

  if(!is.null(clin_file)){


    if(grepl("phenotype",clin_file)){
      clin_xena <- data.table::fread(clin_file,data.table = F)
      index <- intersect(colnames(mrna_expr_xena),clin_xena$submitter_id.samples)
      mrna_xena <- mrna_expr_xena[,index]
      clin_xena <- clin_xena[clin_xena$submitter_id.samples %in% index,]
      clin_xena <- clin_xena[match(colnames(mrna_xena),clin_xena$submitter_id.samples),]
      save(clin_xena, mrna_xena,file = "output_mrna_expr_xena/mrna_and_clin_xena.rdata")

    }else{

      clin_xena <- data.table::fread(clin_file,data.table = F)
      index <- intersect(colnames(mrna_expr_xena),clin_xena$sample)
      mrna_xena <- mrna_expr_xena[,index]
      clin_xena <- clin_xena[clin_xena$sample %in% index,]
      clin_xena <- clin_xena[match(colnames(mrna_xena),clin_xena$sample),]
      save(clin_xena, mrna_xena,file = "output_mrna_expr_xena/mrna_and_clin_xena.rdata")
    }

  }

}
