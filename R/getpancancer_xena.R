#' Prepare pan-cancer expression matrix and clinical information from XENA
#'
#' @description This function can help you prepare the pan-cancer expression
#'    matrix and clinical(phenotype) information.
#' @param tcga_expr_file TCGA pan-cancer gene expression RNA-seq file
#' @param gtex_expr_file GTEx gene expression RNA-seq file
#' @param tcga_clin_file TCGA clin file
#' @param gtex_pheno_file GTEx pheno file
#' @param type type of pan-cancer. "tcga" for only TCGA pan-cancer,
#'    tcga_expr_file and tcga_clin_file should be provided. "gtex" for only GTEx
#'    pan-cancer, gtex_expr_file and gtex_pheno_file should be provided.
#'    "tcga+gtex" for combined TCGA and GTEX pan-cancer, all the 4 files should
#'    be provided.
#' @param save_lnc save lncRNA data, default is "TRUE"
#' @details
#' All the files must be downloaded from the following urls:
#'\itemize{
#'\item{tcga_expr_file: https://toil-xena-hub.s3.us-east-1.amazonaws.com/
#'download/tcga_RSEM_gene_tpm.gz}
#'\item{gtex_expr_file: https://toil-xena-hub.s3.us-east-1.amazonaws.com/
#'download/gtex_RSEM_gene_tpm.gz}
#'\item{tcga_clin_file: https://tcga-pancan-atlas-hub.s3.us-east-1.
#'amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp}
#'\item{gtex_pheno_file: https://toil-xena-hub.s3.us-east-1.amazonaws.com/
#'download/GTEX_phenotype.gz}
#'}
#'
#' @return pan-cancer expression matrix and clinical info. The data are saved
#'    under the directory of "output_pancancer_xena".
#' @importFrom data.table :=
#' @importFrom data.table transpose
#' @importFrom data.table setDT
#' @importFrom data.table setcolorder
#' @export

getpancancer_xena <- function(tcga_expr_file,tcga_clin_file,
                              gtex_expr_file,gtex_pheno_file,
                              type=c("tcga","gtex","tcga+gtex"),
                              save_lnc=TRUE){
  if(!dir.exists("output_pancancer_xena")){dir.create("output_pancancer_xena")}
  if(type=="tcga"){
    message("=> Reading data....")
    tcga_clin <- data.table::fread(tcga_clin_file)
    tcga_expr <- data.table::fread(tcga_expr_file)
    message("=> Cleaning data....")
    index <- c("sample",intersect(tcga_clin$sample, colnames(tcga_expr)[-1]))
    tcga_expr <- tcga_expr[,..index]
    tcga_clin <- tcga_clin[match(colnames(tcga_expr)[-1],tcga_clin$sample),]
    names(tcga_clin)[c(1:3)] <- c("sample_id","patient_id","project")
    names(tcga_expr)[1] <- "gene_id"
    # mRNA
    message("=> Extracting mRNA....")
    tcga_expr_mrna <- merge(gencodev23_mrna,tcga_expr, by = "gene_id")
    tcga_expr_mrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
    tcga_expr_mrna <- tcga_expr_mrna[order(-rmea)]
    tcga_expr_mrna <- tcga_expr_mrna[!duplicated(gene_name)]
    tcga_expr_mrna[,c("gene_id","rmea"):=NULL]
    tcga_expr_mrna <- data.table::transpose(tcga_expr_mrna,make.names = "gene_name",keep.names = NULL)
    tcga_mrna_clin <- cbind(tcga_clin,tcga_expr_mrna)
    if(save_lnc){
      #lncrna
      message("=> Extracting lncRNA....")
      tcga_expr_lncrna <- merge(gencodev23_lncrna,tcga_expr, by = "gene_id")
      tcga_expr_lncrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
      tcga_expr_lncrna <- tcga_expr_lncrna[order(-rmea)]
      tcga_expr_lncrna <- tcga_expr_lncrna[!duplicated(gene_name)]
      tcga_expr_lncrna[,c("gene_id","rmea"):=NULL]
      tcga_expr_lncrna <- data.table::transpose(tcga_expr_lncrna,make.names = "gene_name",keep.names = NULL)
      tcga_lncrna_clin <- cbind(tcga_clin,tcga_expr_lncrna)
      tcga_lncrna_clin <- as.data.frame(tcga_lncrna_clin)
      save(tcga_lncrna_clin,file = "output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
    }
    message("=> Saving data....")
    tcga_mrna_clin <- as.data.frame(tcga_mrna_clin)
    save(tcga_mrna_clin,file = "output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    tcga_expr <- as.data.frame(tcga_expr)
    tcga_clin <- as.data.frame(tcga_clin)
    save(tcga_expr,file = "output_pancancer_xena/TCGA_pancancer_expr.rdata")
    save(tcga_clin,file = "output_pancancer_xena/TCGA_pancancer_clin.rdata")
    message("=> successful. \n")
  }
  else if(type=="gtex"){
    message("=> Reading files....")
    gtex_pheno <- data.table::fread(gtex_pheno_file)
    gtex_expr <- data.table::fread(gtex_expr_file)
    message("=> Cleaning data....")
    names(gtex_pheno)[c(1,3)] <- c("sample_id","primary_site")
    gtex_pheno <- gtex_pheno[!primary_site== "<not provided>",c(1,3)]
    index <- c("sample",intersect(colnames(gtex_expr)[-1],gtex_pheno$sample_id))
    gtex_expr <- gtex_expr[,..index]
    gtex_pheno <- gtex_pheno[match(colnames(gtex_expr)[-1],gtex_pheno$sample_id),]
    names(gtex_expr)[1] <- "gene_id"
    # mRNA
    message("=> Extracting mRNA....")
    gtex_expr_mrna <- merge(gencodev23_mrna,gtex_expr, by = "gene_id")
    gtex_expr_mrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
    gtex_expr_mrna <- gtex_expr_mrna[order(-rmea)]
    gtex_expr_mrna <- gtex_expr_mrna[!duplicated(gene_name)]
    gtex_expr_mrna[,c("gene_id","rmea"):=NULL]
    gtex_expr_mrna <- data.table::transpose(gtex_expr_mrna,make.names = "gene_name",keep.names = NULL)
    gtex_mrna_pheno <- cbind(gtex_pheno,gtex_expr_mrna)
    if(save_lnc){
      # lncRNA
      message("=> Extracting lncRNA....")
      gtex_expr_lncrna <- merge(gencodev23_lncrna,gtex_expr, by = "gene_id")
      gtex_expr_lncrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
      gtex_expr_lncrna <- gtex_expr_lncrna[order(-rmea)]
      gtex_expr_lncrna <- gtex_expr_lncrna[!duplicated(gene_name)]
      gtex_expr_lncrna[,c("gene_id","rmea"):=NULL]
      gtex_expr_lncrna <- data.table::transpose(gtex_expr_lncrna,make.names = "gene_name",keep.names = NULL)
      gtex_lncrna_pheno <- cbind(gtex_pheno,gtex_expr_lncrna)
      gtex_lncrna_pheno <- as.data.frame(gtex_lncrna_pheno)
      save(gtex_lncrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
    }
    message("=> Saving data....")
    gtex_mrna_pheno <- as.data.frame(gtex_mrna_pheno)
    save(gtex_mrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    gtex_pheno <- as.data.frame(gtex_pheno)
    gtex_expr <- as.data.frame(gtex_expr)
    save(gtex_pheno,file = "output_pancancer_xena/GTEx_pancancer_pheno.rdata")
    save(gtex_expr,file = "output_pancancer_xena/GTEx_pancancer_expr.rdata")
    message("=> Successful. \n")
  }
  else {
    message("=> Reading TCGA files....")
    tcga_clin <- data.table::fread(tcga_clin_file)
    tcga_expr <- data.table::fread(tcga_expr_file)
    message("=> Cleaning TCGA data....")
    index <- c("sample",intersect(tcga_clin$sample, colnames(tcga_expr)[-1]))
    tcga_expr <- tcga_expr[,..index]
    tcga_clin <- tcga_clin[match(colnames(tcga_expr)[-1],tcga_clin$sample),]
    names(tcga_clin)[c(1:3)] <- c("sample_id","patient_id","project")
    names(tcga_expr)[1] <- "gene_id"
    # mRNA
    message("=> Extracting TCGA mRNA....")
    tcga_expr_mrna <- merge(gencodev23_mrna,tcga_expr, by = "gene_id")
    tcga_expr_mrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
    tcga_expr_mrna <- tcga_expr_mrna[order(-rmea)]
    tcga_expr_mrna <- tcga_expr_mrna[!duplicated(gene_name)]
    tcga_expr_mrna[,c("gene_id","rmea"):=NULL]
    tcga_expr_mrna <- data.table::transpose(tcga_expr_mrna,make.names = "gene_name",keep.names = NULL)
    tcga_mrna_clin <- cbind(tcga_clin,tcga_expr_mrna)
    if(save_lnc){
      #lncrna
      message("=> Extracting TCGA lncRNA....")
      tcga_expr_lncrna <- merge(gencodev23_lncrna,tcga_expr, by = "gene_id")
      tcga_expr_lncrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
      tcga_expr_lncrna <- tcga_expr_lncrna[order(-rmea)]
      tcga_expr_lncrna <- tcga_expr_lncrna[!duplicated(gene_name)]
      tcga_expr_lncrna[,c("gene_id","rmea"):=NULL]
      tcga_expr_lncrna <- data.table::transpose(tcga_expr_lncrna,make.names = "gene_name",keep.names = NULL)
      tcga_lncrna_clin <- cbind(tcga_clin,tcga_expr_lncrna)
      tcga_lncrna_clin <- as.data.frame(tcga_lncrna_clin)
      save(tcga_lncrna_clin,file = "output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
    }
    message("=> Saving TCGA data....")
    tcga_mrna_clin <- as.data.frame(tcga_mrna_clin)
    save(tcga_mrna_clin,file = "output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    tcga_expr <- as.data.frame(tcga_expr)
    tcga_clin <- as.data.frame(tcga_clin)
    save(tcga_expr,file = "output_pancancer_xena/TCGA_pancancer_expr.rdata")
    save(tcga_clin,file = "output_pancancer_xena/TCGA_pancancer_clin.rdata")
    message("=> Reading GTEx files....")
    gtex_pheno <- data.table::fread(gtex_pheno_file)
    gtex_expr <- data.table::fread(gtex_expr_file)
    message("=> Cleaning GTEx data....")
    names(gtex_pheno)[c(1,3)] <- c("sample_id","primary_site")
    gtex_pheno <- gtex_pheno[!primary_site== "<not provided>",c(1,3)]
    index <- c("sample",intersect(colnames(gtex_expr)[-1],gtex_pheno$sample_id))
    gtex_expr <- gtex_expr[,..index]
    gtex_pheno <- gtex_pheno[match(colnames(gtex_expr)[-1],gtex_pheno$sample_id),]
    names(gtex_expr)[1] <- "gene_id"
    # mRNA
    message("=> Extracting GTEx mRNA....")
    gtex_expr_mrna <- merge(gencodev23_mrna,gtex_expr, by = "gene_id")
    gtex_expr_mrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
    gtex_expr_mrna <- gtex_expr_mrna[order(-rmea)]
    gtex_expr_mrna <- gtex_expr_mrna[!duplicated(gene_name)]
    gtex_expr_mrna[,c("gene_id","rmea"):=NULL]
    gtex_expr_mrna <- data.table::transpose(gtex_expr_mrna,make.names = "gene_name",keep.names = NULL)
    gtex_mrna_pheno <- cbind(gtex_pheno,gtex_expr_mrna)
    if(save_lnc){
      # lncRNA
      message("=> Extracting GTEx lncRNA....")
      gtex_expr_lncrna <- merge(gencodev23_lncrna,gtex_expr, by = "gene_id")
      gtex_expr_lncrna[, `:=`(rmea=rowMeans(.SD)),.SDcols=-c(1:2)]
      gtex_expr_lncrna <- gtex_expr_lncrna[order(-rmea)]
      gtex_expr_lncrna <- gtex_expr_lncrna[!duplicated(gene_name)]
      gtex_expr_lncrna[,c("gene_id","rmea"):=NULL]
      gtex_expr_lncrna <- data.table::transpose(gtex_expr_lncrna,make.names = "gene_name",keep.names = NULL)
      gtex_lncrna_pheno <- cbind(gtex_pheno,gtex_expr_lncrna)
      gtex_lncrna_pheno <- as.data.frame(gtex_lncrna_pheno)
      save(gtex_lncrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
    }
    message("=> Saving GTEx data....")
    gtex_mrna_pheno <- as.data.frame(gtex_mrna_pheno)
    save(gtex_mrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    gtex_pheno <- as.data.frame(gtex_pheno)
    gtex_expr <- as.data.frame(gtex_expr)
    save(gtex_pheno,file = "output_pancancer_xena/GTEx_pancancer_pheno.rdata")
    save(gtex_expr,file = "output_pancancer_xena/GTEx_pancancer_expr.rdata")
    message("=> Combining mRNA....")
    #load("output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    #load("output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    gtex_mrna_pheno <- data.table::setDT(gtex_mrna_pheno)
    tcga_mrna_clin <- data.table::setDT(tcga_mrna_clin)
    index <- gtex_mrna_pheno$primary_site %in% tcga_gtex_match$primary_site
    gtex_mrna_pheno <- gtex_mrna_pheno[index,]
    gtex_mrna_pheno <- merge(gtex_mrna_pheno,tcga_gtex_match,by="primary_site",all.x = F,all.y = F,allow.cartesian=TRUE)
    gtex_mrna_pheno[,`:=`(sample_type="GTEx_normal")]
    tcga_mrna_clin <- tcga_mrna_clin[,c(1,3,35:19759)]
    tcga_mrna_clin[,`:=`(primary_site = project)]
    tcga_mrna_clin[,`:=`(sample_type=ifelse(as.numeric(substr(sample_id,14,15))<10,"TCGA_tumor","TCGA_normal"))]
    data.table::setcolorder(tcga_mrna_clin,match(colnames(tcga_mrna_clin),colnames(gtex_mrna_pheno)))
    tcga_gtex_mrna_pheno <- rbind(tcga_mrna_clin,gtex_mrna_pheno)
    data.table::setcolorder(tcga_gtex_mrna_pheno,c(19728,19729,1,2,3:19727))
    tcga_gtex_mrna_pheno <- as.data.frame(tcga_gtex_mrna_pheno)
    save(tcga_gtex_mrna_pheno,file="output_pancancer_xena/TCGA_GTEx_pancancer_mrna_pheno.rdata")
    if(save_lnc){
      message("=> Combining lncRNA....")
      #load("output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
      #load("output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
      gtex_lncrna_pheno <- data.table::setDT(gtex_lncrna_pheno)
      tcga_lncrna_clin <- data.table::setDT(tcga_lncrna_clin)
      index <- gtex_lncrna_pheno$primary_site %in% tcga_gtex_match$primary_site
      gtex_lncrna_pheno <- gtex_lncrna_pheno[index,]
      gtex_lncrna_pheno <- merge(gtex_lncrna_pheno,tcga_gtex_match,by="primary_site",all.x = F,all.y = F,allow.cartesian=TRUE)
      gtex_lncrna_pheno[,`:=`(sample_type="GTEx_normal")]
      tcga_lncrna_clin <- tcga_lncrna_clin[,c(1,3,35:14395)]
      tcga_lncrna_clin[,`:=`(primary_site = project)]
      tcga_lncrna_clin[,`:=`(sample_type=ifelse(as.numeric(substr(sample_id,14,15))<10,"TCGA_tumor","TCGA_normal"))]
      data.table::setcolorder(tcga_lncrna_clin,match(colnames(tcga_lncrna_clin),colnames(gtex_lncrna_pheno)))
      tcga_gtex_lncrna_pheno <- rbind(tcga_lncrna_clin,gtex_lncrna_pheno)
      data.table::setcolorder(tcga_gtex_lncrna_pheno,c(14364,14365,1,2,3:14363))
      tcga_gtex_lncrna_pheno <- as.data.frame(tcga_gtex_lncrna_pheno)
      save(tcga_gtex_lncrna_pheno,file="output_pancancer_xena/TCGA_GTEx_pancancer_lncrna_pheno.rdata")
    }
    message("=> Successful. \n")
  }

}

