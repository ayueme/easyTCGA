#' Batch survival analysis
#'
#' @description This function can automatically do the univariate cox regression
#'    and log-rank test according to the median of gene expression(no other
#'    choice for current version). Attention: only do survival analysis on tumor
#'    samples
#' @param exprset expression matrix prepared by getmrnaexpr() or getmirnaexpr()
#' @param clin clinical information prepared by getmrnaexpr() or getmirnaexpr()
#' @param expr_type type of expression matrix, one of counts, tpm, fpkm. If
#'     "counts", the expression matrix will be transformed by DESeq2::vst() or
#'     DESeq2::varianceStabilizingTransformation, if "tpm" or "fpkm", it will be
#'     transformed by log2(x + 0.1)
#' @param project a character used as the part of file name
#' @param min_sample_size min sample size of each group for survival analysis,
#'     default is 5
#' @param print_index print index, default is "TRUE"
#'
#' @return a list of survival analysis results, both log-rank and cox regression
#' @export

batch_survival <- function(exprset,
                           clin,
                           expr_type = c("counts","tpm","fpkm"),
                           project,
                           min_sample_size = 5,
                           print_index = TRUE
                           ) {
  cli::cli_alert_info("Only do survival analysis on tumor samples!")
  if (!dir.exists("output_survival")) {
    dir.create("output_survival")
  }

  if(expr_type == "counts"){
    object <- DESeq2::DESeqDataSetFromMatrix(as.matrix(exprset),
                                             S4Vectors::DataFrame(row.names = colnames(exprset)),~1)
    object <- DESeq2::estimateSizeFactors(object)
    nsubs <- sum(rowMeans(DESeq2::counts(object, normalized=TRUE)) > 5 )
    if(nsubs > 1000){
      cli::cli_alert_info("Your exprset will be transformed by DESeq2::vst.")
      exprset <- DESeq2::vst(as.matrix(exprset))
    } else {
      cli::cli_alert_info("Your exprset will be transformed by DESeq2::varianceStabilizingTransformation.")
      exprset <- DESeq2::varianceStabilizingTransformation(as.matrix(exprset))
    }
  } else {
    cli::cli_alert_info("Your exprset will be transformed by log2(x + 0.1).")
    exprset <- log2(exprset + 0.1)
  }
  keep_samples <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
  exprset <- exprset[, keep_samples]
  clin <- clin[keep_samples, c("days_to_last_follow_up", "vital_status")]
  names(clin) <- c("time", "event")
  drop <- is.na(clin$event)
  clin <- clin[!drop, ]
  exprset <- exprset[, !drop]
  clin$event <- ifelse(clin$event == "Dead", 1, 0)

  expr_clin <- cbind(clin, t(exprset))
  gene <- colnames(expr_clin)[-c(1:2)]

  res_survival <- list()

  cox.result <- list()

  for(i in 1:length(gene)){
    if(print_index) print(i)
    group <- ifelse(expr_clin[,gene[i]]>stats::median(expr_clin[,gene[i]]),"high","low")
    if(length(table(group)) == 1) next
    if(length(grep("high",group)) < min_sample_size) next
    if(length(grep("low",group)) < min_sample_size) next
    x <- survival::coxph(survival::Surv(time, event)~group, data = expr_clin)
    tmp1 <- broom::tidy(x,exponentiate = T, conf.int = T)
    cox.result[[i]] <- c(gene[i],tmp1)
  }
  res.cox <- data.frame(do.call(rbind,cox.result))
  names(res.cox)[c(1,3,5)] <- c("gene","HR","Wald_z")

  res_survival[[1]] <- res.cox
  names(res_survival)[[1]] <- "res.cox"

  logrank.result <- list()

  for(i in 1:length(gene)){
    if(print_index) print(i)
    group <- ifelse(expr_clin[,gene[i]]>stats::median(expr_clin[,gene[i]]),"high","low")
    if(length(table(group)) == 1) next
    if(length(grep("high",group)) < min_sample_size) next
    if(length(grep("low",group)) < min_sample_size) next
    x <- survival::survdiff(survival::Surv(time, event)~group, data = expr_clin)
    pValue <- 1-stats::pchisq(x$chisq,df=1)
    logrank.result[[i]] <- c(gene[i],pValue)
  }

  res.logrank <- data.frame(do.call(rbind,logrank.result))
  names(res.logrank) <- c("gene","p.value")

  res_survival[[2]] <- res.logrank
  names(res_survival)[[2]] <- "res.logrank"

  save(res_survival, file = paste0("output_survival/",project,"_survival_results.rdata"))
  cli::cli_alert_success("Analysis done.")
  return(res_survival)

}
