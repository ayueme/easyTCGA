#' Get mRNA and lncRNA expression matrix from XENA file
#'
#' @param file_name gene expression RNAseq(HTSeq-Counts) file downloaded from
#'     https://gdc.xenahubs.net
#'
#' @return mRNA and lncRNA expression matrix. The data are saved in the
#'     directory of "output_mRNA_expr_xena".
#' @export

getmrnaexpr_xena <- function(file_name){

  if(!dir.exists("output_mrna_expr_xena")){
    dir.create("output_mrna_expr_xena")
  }

  exprset <- utils::read.table(file_name,sep = "\t",header = T,check.names = F)
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
  tmp <- 2^tmp-1
  mrna_expr_xena <- floor(tmp)
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
  tmp <- 2^tmp-1
  lncrna_expr_xena <- floor(tmp)
  save(lncrna_expr_xena, file = "output_mrna_expr_xena/lncrna_expr_xena.rdata")
}
