#' Batch survival analysis
#'
#' @param exprset expression set prepared by getmrnaexpr()
#' @param clin clinical information prepared by getmrnaexpr()
#' @param expr_type type of expression set
#' @param project a TCGA project
#' @param min_sample_size min sample size
#' @param print_index print index, default is "TRUE"
#'
#' @return survival analysis resaults
#' @export

globalVariables(c('as.formula','pchisq','na.omit'))
batch_survival <- function(exprset, clin, expr_type=c("counts","tpm","fpkm"),
                           project, min_sample_size=5,
                           print_index = TRUE
) {
  cat("Only do survival analysis on tumor samples! \n")

  if (!file.exists("output_survival")) {
    dir.create("output_survival")
  }

  if(expr_type == "counts"){
    cat("Your exprset will be transformed by DESeq2::vst ")
    exprset <- DESeq2::vst(as.matrix(exprset))
  } else {
    cat("Your exprset will be transformed by log2(x + 0.1) ")
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
    surv <- as.formula(paste('Surv(time, event)~', "group"))
    #tmp <- cbind(expr_clin[,1:2],group)
    x <- survival::coxph(surv, data = expr_clin)
    tmp1 <- broom::tidy(x,exponentiate = T, conf.int = T)
    cox.result[[i]] <- c(gene[i],tmp1)
  }
  res.cox <- data.frame(do.call(rbind,cox.result))
  names(res.cox)[1] <- "gene"

  res_survival[[1]] <- res.cox
  names(res_survival)[[1]] <- "res.cox"

  logrank.result <- list()

  for(i in 1:length(gene)){
    if(print_index) print(i)
    group <- ifelse(expr_clin[,gene[i]]>stats::median(expr_clin[,gene[i]]),"high","low")
    if(length(table(group)) == 1) next
    if(length(grep("high",group)) < min_sample_size) next
    if(length(grep("low",group)) < min_sample_size) next
    surv <- as.formula(paste('Surv(time, event)~', "group"))
    #tmp <- cbind(expr_clin[,1:2],group)
    x <- survival::survdiff(surv, data = expr_clin)
    pValue <- 1-pchisq(x$chisq,df=1)
    logrank.result[[i]] <- c(gene[i],pValue)
  }

  res.logrank <- data.frame(do.call(rbind,logrank.result))
  names(res.logrank) <- c("gene","p.value")

  res_survival[[2]] <- res.logrank
  names(res_survival)[[2]] <- "res.logrank"

  save(res_survival, file = paste0("output_survival/",project,"_survival_results.rdata"))
  return(res_survival)
}
