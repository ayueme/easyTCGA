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
#' @export

getpancancer_xena <- function(tcga_expr_file,tcga_clin_file,
                              gtex_expr_file,gtex_pheno_file,
                              type=c("tcga","gtex","tcga+gtex")){

  if(!dir.exists("output_pancancer_xena")){dir.create("output_pancancer_xena")}

  if(type=="tcga"){

    message("=> Reading files....")
    # clean data
    tcga_clin <- data.table::fread(tcga_clin_file,data.table = F)
    tcga_expr <- data.table::fread(tcga_expr_file,data.table = F)
    message("=> Cleaning data....")
    index <- intersect(tcga_clin$sample, colnames(tcga_expr)[-1])
    tcga_expr <- tcga_expr[,c("sample",index)]
    tcga_clin <- tcga_clin[match(colnames(tcga_expr)[-1],tcga_clin$sample),]
    names(tcga_clin)[c(1:3)] <- c("sample_id","patient_id","project")
    names(tcga_expr)[1] <- "gene_id"
    save(tcga_expr,file = "output_pancancer_xena/TCGA_pancancer_expr.rdata")
    save(tcga_clin,file = "output_pancancer_xena/TCGA_pancancer_clin.rdata")
    # mRNA
    message("=> Extracting mRNA....")
    tcga_expr_mrna <- merge(gencodev23_mrna,tcga_expr, by = "gene_id")
    tcga_expr_mrna <- tcga_expr_mrna[,-1]
    rowm <- as.data.frame(rowMeans(tcga_expr_mrna[, -1]))
    names(rowm) <- "rmea"
    tcga_expr_mrna <- cbind(rowm, tcga_expr_mrna)
    tcga_expr_mrna <- tcga_expr_mrna[order(tcga_expr_mrna$rmea,decreasing=T),]
    tcga_expr_mrna <- tcga_expr_mrna[!duplicated(tcga_expr_mrna$gene_name),]
    rownames(tcga_expr_mrna) <- tcga_expr_mrna[,2]
    tcga_expr_mrna <- tcga_expr_mrna[,c(-1,-2)]
    tcga_mrna_clin <- cbind(tcga_clin,as.data.frame(t(tcga_expr_mrna)))
    save(tcga_mrna_clin,file = "output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    #lncrna
    rm(tcga_expr_mrna);rm(tcga_mrna_clin)
    message("=> Extracting lncRNA....")
    tcga_expr_lncrna <- merge(gencodev23_lncrna,tcga_expr, by = "gene_id")
    tcga_expr_lncrna <- tcga_expr_lncrna[,-1]
    rowm <- as.data.frame(rowMeans(tcga_expr_lncrna[, -1]))
    names(rowm) <- "rmea"
    tcga_expr_lncrna <- cbind(rowm, tcga_expr_lncrna)
    tcga_expr_lncrna <- tcga_expr_lncrna[order(tcga_expr_lncrna$rmea,decreasing=T),]
    tcga_expr_lncrna <- tcga_expr_lncrna[!duplicated(tcga_expr_lncrna$gene_name),]
    rownames(tcga_expr_lncrna) <- tcga_expr_lncrna[,2]
    tcga_expr_lncrna <- tcga_expr_lncrna[,c(-1,-2)]
    tcga_lncrna_clin <- cbind(tcga_clin,as.data.frame(t(tcga_expr_lncrna)))
    save(tcga_lncrna_clin,file = "output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
    message("=> Successful. \n")
  }
  else if(type=="gtex"){

    message("=> Reading files....")
    ## clean data
    gtex_pheno <- utils::read.table(gtex_pheno_file,sep = "\t",header = T)
    gtex_pheno <- subset(gtex_pheno, X_primary_site != "<not provided>",select = c(1,3))
    names(gtex_pheno) <- c("sample_id","primary_site")
    gtex_expr <- data.table::fread(gtex_expr_file,data.table = F)
    message("=> Cleaning data....")
    index <- intersect(colnames(gtex_expr)[-1],gtex_pheno$sample_id)
    gtex_expr <- gtex_expr[,c("sample",index)]
    gtex_pheno <- gtex_pheno[match(colnames(gtex_expr)[-1],gtex_pheno$sample_id),]
    names(gtex_expr)[1] <- "gene_id"
    save(gtex_pheno,file = "output_pancancer_xena/GTEx_pancancer_pheno.rdata")
    save(gtex_expr,file = "output_pancancer_xena/GTEx_pancancer_expr.rdata")
    # mRNA
    message("=> Extracting mRNA....")
    gtex_expr_mrna <- merge(gencodev23_mrna,gtex_expr, by = "gene_id")
    gtex_expr_mrna <- gtex_expr_mrna[,-1]
    rowm <- as.data.frame(rowMeans(gtex_expr_mrna[, -1]))
    names(rowm) <- "rmea"
    gtex_expr_mrna <- cbind(rowm, gtex_expr_mrna)
    gtex_expr_mrna <- gtex_expr_mrna[order(gtex_expr_mrna$rmea,decreasing=T),]
    gtex_expr_mrna <- gtex_expr_mrna[!duplicated(gtex_expr_mrna$gene_name),]
    rownames(gtex_expr_mrna) <- gtex_expr_mrna[,2]
    gtex_expr_mrna <- gtex_expr_mrna[,c(-1,-2)]
    gtex_mrna_pheno <- cbind(gtex_pheno,as.data.frame(t(gtex_expr_mrna)))
    save(gtex_mrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    # lncRNA
    rm(gtex_mrna_pheno);rm(gtex_expr_mrna)
    message("=> Extracting lncRNA....")
    gtex_expr_lncrna <- merge(gencodev23_lncrna,gtex_expr, by = "gene_id")
    gtex_expr_lncrna <- gtex_expr_lncrna[,-1]
    rowm <- as.data.frame(rowMeans(gtex_expr_lncrna[, -1]))
    names(rowm) <- "rmea"
    gtex_expr_lncrna <- cbind(rowm, gtex_expr_lncrna)
    gtex_expr_lncrna <- gtex_expr_lncrna[order(gtex_expr_lncrna$rmea,decreasing=T),]
    gtex_expr_lncrna <- gtex_expr_lncrna[!duplicated(gtex_expr_lncrna$gene_name),]
    rownames(gtex_expr_lncrna) <- gtex_expr_lncrna[,2]
    gtex_expr_lncrna <- gtex_expr_lncrna[,c(-1,-2)]
    gtex_lncrna_pheno <- cbind(gtex_pheno,as.data.frame(t(gtex_expr_lncrna)))
    save(gtex_lncrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
    message("=> Successful. \n")
  }
  else {

    message("=> Preparing TCGA data....")
    ### clean data
    ## TCGA
    message(" => Reading data....")
    tcga_clin <- data.table::fread(tcga_clin_file,data.table = F)
    tcga_expr <- data.table::fread(tcga_expr_file,data.table = F)
    message(" => Cleaning data....")
    index <- intersect(tcga_clin$sample, colnames(tcga_expr)[-1])
    tcga_expr <- tcga_expr[,c("sample",index)]
    tcga_clin <- tcga_clin[match(colnames(tcga_expr)[-1],tcga_clin$sample),]
    names(tcga_clin)[c(1:3)] <- c("sample_id","patient_id","project")
    names(tcga_expr)[1] <- "gene_id"
    save(tcga_expr,file = "output_pancancer_xena/TCGA_pancancer_expr.rdata")
    save(tcga_clin,file = "output_pancancer_xena/TCGA_pancancer_clin.rdata")
    # mRNA
    message(" => Extracting mRNA....")
    tcga_expr_mrna <- merge(gencodev23_mrna,tcga_expr, by = "gene_id")
    tcga_expr_mrna <- tcga_expr_mrna[,-1]
    rowm <- as.data.frame(rowMeans(tcga_expr_mrna[, -1]))
    names(rowm) <- "rmea"
    tcga_expr_mrna <- cbind(rowm, tcga_expr_mrna)
    tcga_expr_mrna <- tcga_expr_mrna[order(tcga_expr_mrna$rmea,decreasing=T),]
    tcga_expr_mrna <- tcga_expr_mrna[!duplicated(tcga_expr_mrna$gene_name),]
    rownames(tcga_expr_mrna) <- tcga_expr_mrna[,2]
    tcga_expr_mrna <- tcga_expr_mrna[,c(-1,-2)]
    tcga_mrna_clin <- cbind(tcga_clin,as.data.frame(t(tcga_expr_mrna)))
    save(tcga_mrna_clin,file = "output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    #lncrna
    rm(tcga_expr_mrna);rm(tcga_mrna_clin)
    message(" => Extracting lncRNA....")
    tcga_expr_lncrna <- merge(gencodev23_lncrna,tcga_expr, by = "gene_id")
    tcga_expr_lncrna <- tcga_expr_lncrna[,-1]
    rowm <- as.data.frame(rowMeans(tcga_expr_lncrna[, -1]))
    names(rowm) <- "rmea"
    tcga_expr_lncrna <- cbind(rowm, tcga_expr_lncrna)
    tcga_expr_lncrna <- tcga_expr_lncrna[order(tcga_expr_lncrna$rmea,decreasing=T),]
    tcga_expr_lncrna <- tcga_expr_lncrna[!duplicated(tcga_expr_lncrna$gene_name),]
    rownames(tcga_expr_lncrna) <- tcga_expr_lncrna[,2]
    tcga_expr_lncrna <- tcga_expr_lncrna[,c(-1,-2)]
    tcga_lncrna_clin <- cbind(tcga_clin,as.data.frame(t(tcga_expr_lncrna)))
    save(tcga_lncrna_clin,file = "output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
    #save(tcga_expr_lncrna,file = "output_pancancer_xena/TCGA_pancancer_expr_lncrna.rdata")

    ## gteX
    message("=> Preparing GTEx data....")
    rm(list = ls())
    message(" => Reading data....")
    gtex_pheno <- utils::read.table(gtex_pheno_file,sep = "\t",header = T)
    gtex_pheno <- subset(gtex_pheno, X_primary_site != "<not provided>",select = c(1,3))
    names(gtex_pheno) <- c("sample_id","primary_site")
    gtex_expr <- data.table::fread(gtex_expr_file,data.table = F)
    message(" => Cleaning data....")
    index <- intersect(colnames(gtex_expr)[-1],gtex_pheno$sample_id)
    gtex_expr <- gtex_expr[,c("sample",index)]
    gtex_pheno <- gtex_pheno[match(colnames(gtex_expr)[-1],gtex_pheno$sample_id),]
    names(gtex_expr)[1] <- "gene_id"
    save(gtex_pheno,file = "output_pancancer_xena/GTEx_pancancer_pheno.rdata")
    save(gtex_expr,file = "output_pancancer_xena/GTEx_pancancer_expr.rdata")
    # mRNA
    message(" => Extracting mRNA....")
    gtex_expr_mrna <- merge(gencodev23_mrna,gtex_expr, by = "gene_id")
    gtex_expr_mrna <- gtex_expr_mrna[,-1]
    rowm <- as.data.frame(rowMeans(gtex_expr_mrna[, -1]))
    names(rowm) <- "rmea"
    gtex_expr_mrna <- cbind(rowm, gtex_expr_mrna)
    gtex_expr_mrna <- gtex_expr_mrna[order(gtex_expr_mrna$rmea,decreasing=T),]
    gtex_expr_mrna <- gtex_expr_mrna[!duplicated(gtex_expr_mrna$gene_name),]
    rownames(gtex_expr_mrna) <- gtex_expr_mrna[,2]
    gtex_expr_mrna <- gtex_expr_mrna[,c(-1,-2)]
    gtex_mrna_pheno <- cbind(gtex_pheno,as.data.frame(t(gtex_expr_mrna)))
    save(gtex_mrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    # lncRNA
    rm(gtex_mrna_pheno);rm(gtex_expr_mrna)
    message(" => Extracting lncRNA....")
    gtex_expr_lncrna <- merge(gencodev23_lncrna,gtex_expr, by = "gene_id")
    gtex_expr_lncrna <- gtex_expr_lncrna[,-1]
    rowm <- as.data.frame(rowMeans(gtex_expr_lncrna[, -1]))
    names(rowm) <- "rmea"
    gtex_expr_lncrna <- cbind(rowm, gtex_expr_lncrna)
    gtex_expr_lncrna <- gtex_expr_lncrna[order(gtex_expr_lncrna$rmea,decreasing=T),]
    gtex_expr_lncrna <- gtex_expr_lncrna[!duplicated(gtex_expr_lncrna$gene_name),]
    rownames(gtex_expr_lncrna) <- gtex_expr_lncrna[,2]
    gtex_expr_lncrna <- gtex_expr_lncrna[,c(-1,-2)]
    gtex_lncrna_pheno <- cbind(gtex_pheno,as.data.frame(t(gtex_expr_lncrna)))
    save(gtex_lncrna_pheno,file = "output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
    # 整合TCGA 和 GTEx
    message("=> Combining TCGA and GTEx....")
    rm(list = ls())
    # 整合mRNA
    message(" => Combining mRNA....")
    load("output_pancancer_xena/TCGA_pancancer_mrna_clin.rdata")
    load("output_pancancer_xena/GTEx_pancancer_mrna_pheno.rdata")
    #只选在tcga里有对应的样本
    index <- gtex_mrna_pheno$primary_site %in% tcga_gtex_match$primary_site
    gtex_mrna_pheno <- gtex_mrna_pheno[index,]
    gtex_mrna_pheno <- merge(gtex_mrna_pheno,tcga_gtex_match,by="primary_site",all.x = F,all.y = F)
    gtex_mrna_pheno$sample_type <- "GTEx_normal"
    # tcga里面前34列都是临床信息，我们只要2列即可
    tcga_mrna_clin <- tcga_mrna_clin[,c(1,3,35:ncol(tcga_mrna_clin))]
    tcga_mrna_clin$primary_site <- tcga_mrna_clin$project
    tcga_mrna_clin$sample_type <- ifelse(as.numeric(substr(tcga_mrna_clin$sample_id,14,15))<10,"TCGA_tumor","TCGA_normal")
    gtex_mrna_pheno <- gtex_mrna_pheno[,match(colnames(tcga_mrna_clin),colnames(gtex_mrna_pheno))]
    tcga_gtex_mrna_pheno <- rbind(tcga_mrna_clin,gtex_mrna_pheno)
    tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[,c(1,2,19728,19729,3:19727)]
    save(tcga_gtex_mrna_pheno,file="output_pancancer_xena/TCGA_GTEx_pancancer_mrna_pheno.rdata")
    # 整合lncRNA
    rm(list = ls())
    message(" => Combining lncRNA....")
    load("output_pancancer_xena/GTEx_pancancer_lncrna_pheno.rdata")
    load("output_pancancer_xena/TCGA_pancancer_lncrna_clin.rdata")
    index <- gtex_lncrna_pheno$primary_site %in% tcga_gtex_match$primary_site
    gtex_lncrna_pheno <- gtex_lncrna_pheno[index,]
    gtex_lncrna_pheno <- merge(gtex_lncrna_pheno,tcga_gtex_match,by="primary_site",all.x = F,all.y = F)
    gtex_lncrna_pheno$sample_type <- "GTEx_normal"
    tcga_lncrna_clin <- tcga_lncrna_clin[,c(1,3,35:ncol(tcga_lncrna_clin))]
    tcga_lncrna_clin$primary_site <- tcga_lncrna_clin$project
    tcga_lncrna_clin$sample_type <- ifelse(as.numeric(substr(tcga_lncrna_clin$sample_id,14,15))<10,"TCGA_tumor","TCGA_normal")
    gtex_lncrna_pheno <- gtex_lncrna_pheno[,match(colnames(tcga_lncrna_clin),colnames(gtex_lncrna_pheno))]
    tcga_gtex_lncrna_pheno <- rbind(tcga_lncrna_clin,gtex_lncrna_pheno)
    tcga_gtex_lncrna_pheno <- tcga_gtex_lncrna_pheno[,c(1,2,14364,14365,3:14363)]
    save(tcga_gtex_lncrna_pheno,file="output_pancancer_xena/TCGA_GTEx_pancancer_lncrna_pheno.rdata")
    message("=> successful. \n")
  }


}
