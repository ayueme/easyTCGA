#' Box plot of marker expression between groups
#'
#' @param exprset expression matrix
#' @param marker marker(s) you want to display
#' @param group the groups of your samples
#'
#' @return box plot and the plotting data
#' @export
#'

plot_gene <- function(exprset,marker,group){
  tmp <- t(exprset)
  tmp <- as.data.frame(tmp[,marker])
  tmp$sample_id <- rownames(tmp)

  # only 1 marker
  if(length(marker)<2){
    names(tmp)[1] <- "expression"
    p1 <- ggplot2::ggplot(tmp,ggplot2::aes(group,expression,fill=group))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x=NULL,y=marker)+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "none")+
      ggpubr::stat_compare_means(ggplot2::aes(group = group,label = "p.format"),
                                 method = "kruskal.test")
    print(p1)
    return(tmp)
  } else {
    # more than 1 marker
    tmp$group <- group
    tmp <- reshape2::melt(tmp, id.vars=c("sample_id","group"),
                          variable.name="markers",value.name = "expression")
    p1 <- ggplot2::ggplot(tmp,ggplot2::aes(group,expression,fill=group))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(x=NULL,y="expression")+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "none"
                     ,axis.text.x = ggplot2::element_text(angle = 45,hjust = 1))+
      ggpubr::stat_compare_means(ggplot2::aes(group = group,label = "p.format"),
                                 method = "kruskal.test")+
      ggplot2::facet_wrap(~markers, scales = "free_y")
    print(p1)
    return(tmp)
  }
}




#' Box plot of marker expression between paired normal and tumor samples
#'
#' @param exprset expression matrix
#' @param marker marker(s) you want to display
#'
#' @return box plot and the plotting data
#' @export
#'
globalVariables("sample_id")
plot_gene_paired <- function(exprset, marker){
  # get paired samples
  sample_group <- ifelse(as.numeric(substr(colnames(exprset),14,15))<10,"tumor","normal")
  tmp <- data.frame(sample_group = sample_group, sample_id=colnames(exprset))
  tmp_nor <- tmp[tmp$sample_group=="normal",]
  tmp_tum <- tmp[tmp$sample_group=="tumor",]
  patient <- substr(tmp_nor$sample_id,1,12)
  keep <- substr(tmp_tum$sample_id,1,12) %in% patient
  tmp_tum <- tmp_tum[keep,]
  tmp_tum <- tmp_tum[!duplicated(substr(tmp_tum$sample_id,1,12)),]
  tmp_pair <- rbind(tmp_tum,tmp_nor)
  # get marker expr
  tmp <- t(exprset)
  tmp <- as.data.frame(tmp[,marker])
  tmp$sample_id <- rownames(tmp)

  if(length(marker)<2){
    # only one marker
    names(tmp)[1] <- "expression"
    tmp_pair <- merge(tmp_pair, tmp, by="sample_id")
    tmp_pair$sample_id <- substr(tmp_pair$sample_id,1,12)

    p2 <- ggplot2::ggplot(tmp_pair,ggplot2::aes(sample_group,expression,color=sample_group))+
      ggplot2::geom_boxplot()+
      ggplot2::geom_point(size=3)+
      ggplot2::geom_line(ggplot2::aes(group=sample_id),color="grey")+
      ggplot2::scale_color_manual(values = c("#028EA1","#F2AA9D"))+
      ggplot2::labs(x=NULL,y=paste0("paired_",marker))+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "none")+
      ggpubr::stat_compare_means(ggplot2::aes(group = sample_group,label = "p.format"),
                                 method = "kruskal.test",paired = T)
    print(p2)
    return(tmp_pair)

  } else {
    # more than one marker
    #names(tmp) <- marker
    tmp_pair <- merge(tmp_pair, tmp, by="sample_id")
    tmp_pair$sample_id <- substr(tmp_pair$sample_id,1,12)

    tmp_pair <- reshape2::melt(tmp_pair, id.vars=c("sample_id","sample_group"),
                               variable.name="markers",value.name = "expression")
    p2 <- ggplot2::ggplot(tmp_pair,ggplot2::aes(sample_group,expression,color=sample_group))+
      ggplot2::geom_boxplot()+
      ggplot2::geom_point(size=3)+
      ggplot2::geom_line(ggplot2::aes(group=sample_id),color="grey")+
      ggplot2::scale_color_manual(values = c("#028EA1","#F2AA9D"))+
      ggplot2::labs(x=NULL,y="paired_expression")+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "none")+
      ggpubr::stat_compare_means(ggplot2::aes(group = sample_group,label = "p.format"),
                                 method = "kruskal.test",paired = T)+
      ggplot2::facet_wrap(~markers, scales = "free_y")

    print(p2)
    return(tmp_pair)

  }
}







#' K-M plot according to the expression of marker
#'
#' @param exprset expression matrix
#' @param marker marker you want to display
#' @param clin a data.frame with two columns: "time" and "event", and 1 for
#'    live, 0 for dead.
#' @param optimal_cut use "optimal" cutpoint to do survival analysis. If FALSE,
#'    median of expression will be used. Optimal cutpoint is calculated by
#'    survminer::surv_cutpoint().
#'
#' @return K-M plot and plotting data
#' @export
#'

plot_KM <- function(exprset, marker, clin,optimal_cut=TRUE){
  exprset <- exprset[marker,]
  keep_samples <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
  exprset <- exprset[, keep_samples]
  clin <- clin[keep_samples, ]
  #names(clin) <- c("time", "event")
  drop <- is.na(clin$event)
  clin <- clin[!drop, ]
  exprset <- exprset[, !drop]
  #clin$event <- ifelse(clin$event == "Dead", 1, 0)
  expr_clin <- cbind(clin, t(exprset))
  colnames(expr_clin)[3] <- "gene"

  if(optimal_cut){
    #exprset <- exprset[apply(exprset,1,function(x) sum(duplicated(x)) < 0.8*ncol(exprset)),]
    res.cut <- survminer::surv_cutpoint(expr_clin, time = "time", event = "event",
                                        variables = "gene", minprop = 0.0001,progressbar = F)
    res.cat <- survminer::surv_categorize(res.cut)
    fit <- survival::survfit(survival::Surv(time, event) ~ gene, data = res.cat)
    p1 <- survminer::ggsurvplot(fit, data = res.cat, pval = T)
    print(p1)
    res <- list(surv_df = expr_clin, cut=res.cut[["cutpoint"]][1])
    return(res)
  } else {
    expr_clin$group <- ifelse(expr_clin$gene > stats::median(expr_clin$gene),"high","low")
    fit <- survival::survfit(survival::Surv(time, event) ~ group, data = expr_clin)
    p1 <- survminer::ggsurvplot(fit, data = expr_clin, pval = T)
    print(p1)
    return(expr_clin)
  }
}


