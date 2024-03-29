#' Batch survival analysis
#'
#' @description This function can automatically do log-rank test and univariate
#' cox regression according to the "optimal" cutpoint or median of gene expression.
#' Support mRNA, lncRNA and miRNA expression matrix(count, tpm, fpkm, rpm are all
#' available).
#' @param exprset expression matrix from [getmrnaexpr()] or [getmirnaexpr()].
#' @param clin clinical data from [getmrnaexpr()] or [getmirnaexpr()].
#' @param is_count is the expression matrix count? Default is FALSE. If TRUE,
#' the expression matrix will be transformed by [DESeq2::vst()] or
#' [DESeq2::varianceStabilizingTransformation()]. For other types, the
#' expression matrix will be detected automatically to decide whether to be
#' transformed by log2(x + 0.1).
#' @param optimal_cut use "optimal" cutpoint to do survival analysis? Default is
#' TRUE.If FALSE, median of gene expression will be used. Optimal cutpoint is
#' calculated by [survminer::surv_cutpoint()].
#' @param project characters used as the part of file name when saving results.
#' @param save save results to local, default is FALSE.
#' @param min_sample_size min sample size of each group for survival analysis,
#' default is 5.
#' @param print_index print gene index, default is FALSE.
#'
#' @return a list of survival analysis results, including both log-rank and cox
#' regression.
#' @export

batch_survival <- function(exprset, clin, is_count = FALSE, optimal_cut = TRUE,
                           project = NULL, save = FALSE, min_sample_size = 5,
                           print_index = FALSE) {
  # 去除没有生存结局的，需要样本完全对应
  clin <- clin[, c("days_to_last_follow_up", "days_to_death","vital_status")]
  clin$time <- ifelse(clin$vital_status == "Dead", clin$days_to_death,clin$days_to_last_follow_up)
  clin$event = ifelse(clin$vital_status == "Dead",1,0)
  clin <- clin[,c("time","event")]
  drop <- is.na(clin$event)
  clin <- clin[!drop, ]
  exprset <- exprset[, !drop]
  message(paste0("=> Removed ",sum(drop)," samples without events."))

  #去除normal样本
  drop_normal <- as.numeric(substr(colnames(exprset), 14, 15))>9
  exprset <- exprset[, !drop_normal]
  clin <- clin[!drop_normal,]
  message(paste0("=> Removed ",sum(drop_normal)," normal samples."))

  # 表达矩阵处理
  if (is_count) {
    object <- DESeq2::DESeqDataSetFromMatrix(
      as.matrix(exprset), S4Vectors::DataFrame(row.names = colnames(exprset)), ~1)
    object <- DESeq2::estimateSizeFactors(object)
    nsubs <- sum(rowMeans(DESeq2::counts(object, normalized = TRUE)) > 5)
    if (nsubs > 1000) {
      exprset <- DESeq2::vst(as.matrix(exprset))
      message("=> exprset transformed by DESeq2::vst")
    } else {
      message("=> exprset transformed by DESeq2::varianceStabilizingTransformation")
      exprset <- DESeq2::varianceStabilizingTransformation(as.matrix(exprset),fitType = "local")
    }
  } else {
    exprset <- limma::normalizeBetweenArrays(exprset)

    ex <- exprset
    qx <- as.numeric(stats::quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

    if (LogC) {
      # ex[which(ex <= 0)] <- NaN
      exprset <- log2(ex + 0.1)
      message("=> exprset transformed by log2(x + 0.1)")
    }#else{message("=> log2 transform not needed")}
  }

  # 批量生存分析
  res_survival <- list()
  logrank.result <- list()
  cox.result <- list()

  # 最佳截点的生存分析
  if (optimal_cut) {
    # 表达矩阵处理，基因名有-,在计算最佳截点时报错，-会变为.，所以这里先改一下
    aa <- length(unlist(strsplit(rownames(exprset)[1],"-")))>2
    if(aa){rownames(exprset) <- gsub("-",".",rownames(exprset))}

    message("=> Finding optimal cutpoint...")
    exprset <- exprset[apply(exprset,1,function(x)sum(duplicated(x))<0.8*ncol(exprset)),]
    expr_clin <- cbind(clin, t(exprset))
    gene <- colnames(expr_clin)[-c(1:2)]
    #mirna由于名字问题导致报错，-会变为.
    res.cut <- survminer::surv_cutpoint(expr_clin, time = "time", event = "event",
                                        variables = gene, minprop = 0.0001,progressbar = F)
    res.cat <- survminer::surv_categorize(res.cut)

    message("=> Running batch log-rank...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- res.cat[, i + 2]
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::survdiff(survival::Surv(time, event) ~ group, data = res.cat)
      pValue <- 1 - stats::pchisq(x$chisq, df = 1)
      logrank.result[[i]] <- c(gene[i], pValue, res.cut[["cutpoint"]][i, 1])
    }
    res.logrank <- data.frame(do.call(rbind, logrank.result))
    names(res.logrank) <- c("gene", "p.value", "cutpoint")
    res.logrank$p.value <- as.numeric(res.logrank$p.value)
    res.logrank$cutpoint <- as.numeric(res.logrank$cutpoint)
    # 这里再把基因名改回来：
    if(aa){res.logrank$gene <- gsub("\\.","-",res.logrank$gene)}

    res_survival[[1]] <- res.logrank
    names(res_survival)[[1]] <- "res.logrank"

    message("=> Running batch Cox...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- res.cat[, i + 2]
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::coxph(survival::Surv(time, event) ~ group, data = res.cat)
      tmp1 <- broom::tidy(x, exponentiate = T, conf.int = T)
      cox.result[[i]] <- cbind(gene[i], tmp1, res.cut[["cutpoint"]][i, 1])
    }
    res.cox <- data.frame(do.call(rbind, cox.result))
    names(res.cox)[c(1, 3, 5, 9)] <- c("gene", "HR", "Wald_z", "cutpoint")
    # 这里再把基因名改回来：
    if(aa){res.cox$gene <- gsub("\\.","-",res.cox$gene)}

    res_survival[[2]] <- res.cox
    names(res_survival)[[2]] <- "res.cox"

    if(save){
      if (!dir.exists("output_survival")){dir.create("output_survival")}
      if(is.null(project))stop("project should be provided!")
      save(res_survival,file=paste0("output_survival/",project,"_survival_results.rdata"))
    }
    message("=> Analysis done.")
    return(res_survival)
  }
  else { # 不需要最佳截点的生存分析
    expr_clin <- cbind(clin, t(exprset))
    gene <- colnames(expr_clin)[-c(1:2)]

    message("=> Running batch log-rank...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- ifelse(expr_clin[, gene[i]] > stats::median(expr_clin[, gene[i]]), "high", "low")
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::survdiff(survival::Surv(time, event) ~ group, data = expr_clin)
      pValue <- 1 - stats::pchisq(x$chisq, df = 1)
      logrank.result[[i]] <- c(gene[i], pValue)
    }
    res.logrank <- data.frame(do.call(rbind, logrank.result))
    names(res.logrank) <- c("gene", "p.value")
    res.logrank$p.value <- as.numeric(res.logrank$p.value)

    res_survival[[1]] <- res.logrank
    names(res_survival)[[1]] <- "res.logrank"

    message("=> Running batch Cox...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- ifelse(expr_clin[, gene[i]] > stats::median(expr_clin[, gene[i]]), "high", "low")
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::coxph(survival::Surv(time, event) ~ group, data = expr_clin)
      tmp1 <- broom::tidy(x, exponentiate = T, conf.int = T)
      cox.result[[i]] <- cbind(gene[i], tmp1)
    }
    res.cox <- data.frame(do.call(rbind, cox.result))
    names(res.cox)[c(1, 3, 5)] <- c("gene", "HR", "Wald_z")

    res_survival[[2]] <- res.cox
    names(res_survival)[[2]] <- "res.cox"

    if(save){
      if (!dir.exists("output_survival")){dir.create("output_survival")}
      if(is.null(project))stop("project should be provided!")
      save(res_survival,file=paste0("output_survival/",project,"_survival_results.rdata"))
    }
    message("=> Analysis done.")
    return(res_survival)
  }
}
