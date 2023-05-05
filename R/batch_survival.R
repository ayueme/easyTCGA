#' Batch survival analysis
#'
#' @description This function can automatically do log-rank test and univariate
#'    cox regression according to the "optimal" cutpoint or median of gene
#'    expression. Support mRNA, lncRNA and miRNA expression matrix(count, tpm,
#'    fpkm are all supported).
#' @param exprset Expression matrix from getmrnaexpr() or getmirnaexpr().
#' @param clin Clinical information from getmrnaexpr() or getmirnaexpr().
#' @param expr_type Type of expression matrix, one of count, tpm, fpkm.
#'    If "count", the expression matrix will be transformed by DESeq2::vst() or
#'    DESeq2::varianceStabilizingTransformation. For other types, the expression
#'    will be transformed by log2(x + 0.1).
#' @param optimal_cut Use "optimal" cutpoint to do survival analysis. If FALSE,
#'    median of gene expression will be used. Optimal cutpoint is calculated by
#'    survminer::surv_cutpoint().
#' @param project Characters used as the part of file name when saving results.
#' @param min_sample_size Min sample size of each group for survival analysis,
#'    default is 5.
#' @param print_index Print index, default is TRUE.
#'
#' @return a list of survival analysis results, both log-rank and cox regression
#' @export

batch_survival <- function(exprset,
                           clin,
                           expr_type = c("counts", "tpm", "fpkm"),
                           optimal_cut = TRUE,
                           project,
                           min_sample_size = 5,
                           print_index = TRUE) {
  if (!dir.exists("output_survival")) {
    dir.create("output_survival")
  }

  if (expr_type == "counts") {
    object <- DESeq2::DESeqDataSetFromMatrix(
      as.matrix(exprset),
      S4Vectors::DataFrame(row.names = colnames(exprset)), ~1
    )
    object <- DESeq2::estimateSizeFactors(object)
    nsubs <- sum(rowMeans(DESeq2::counts(object, normalized = TRUE)) > 5)
    if (nsubs > 1000) {
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
  cli::cli_alert_info("Only do survival analysis on tumor samples!")
  keep_samples <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
  exprset <- exprset[, keep_samples]
  clin <- clin[keep_samples, c("days_to_last_follow_up", "vital_status")]
  names(clin) <- c("time", "event")
  drop <- is.na(clin$event)
  clin <- clin[!drop, ]
  exprset <- exprset[, !drop]
  clin$event <- ifelse(clin$event == "Dead", 1, 0)

  res_survival <- list()
  logrank.result <- list()
  cox.result <- list()

  if (optimal_cut) {
    cli::cli_alert_info(" => Finding optimal cutpoint.")
    exprset <- exprset[apply(exprset, 1, function(x) sum(duplicated(x)) < 0.8 * ncol(exprset)), ]
    expr_clin <- cbind(clin, t(exprset))
    gene <- colnames(expr_clin)[-c(1:2)]

    res.cut <- survminer::surv_cutpoint(expr_clin,
      time = "time", event = "event",
      variables = gene, minprop = 0.0001, progressbar = F
    )
    res.cat <- survminer::surv_categorize(res.cut)

    cli::cli_alert_info(" => Running batch log-rank...")
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

    res_survival[[1]] <- res.logrank
    names(res_survival)[[1]] <- "res.logrank"

    cli::cli_alert_info(" => Running batch Cox...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- res.cat[, i + 2]
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::coxph(survival::Surv(time, event) ~ group, data = res.cat)
      tmp1 <- broom::tidy(x, exponentiate = T, conf.int = T)
      cox.result[[i]] <- c(gene[i], tmp1, res.cut[["cutpoint"]][i, 1])
    }
    res.cox <- data.frame(do.call(rbind, cox.result))
    names(res.cox)[c(1, 3, 5, 9)] <- c("gene", "HR", "Wald_z", "cutpoint")

    res_survival[[2]] <- res.cox
    names(res_survival)[[2]] <- "res.cox"

    save(res_survival, file = paste0("output_survival/", project, "_survival_results.rdata"))
    cli::cli_alert_success("Analysis done.")
    return(res_survival)
  } else {
    expr_clin <- cbind(clin, t(exprset))
    gene <- colnames(expr_clin)[-c(1:2)]

    cli::cli_alert_info(" => Running batch log-rank...")
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

    res_survival[[1]] <- res.logrank
    names(res_survival)[[1]] <- "res.logrank"

    cli::cli_alert_info(" => Running batch Cox...")
    for (i in 1:length(gene)) {
      if (print_index) print(i)
      group <- ifelse(expr_clin[, gene[i]] > stats::median(expr_clin[, gene[i]]), "high", "low")
      if (length(table(group)) == 1) next
      if (length(grep("high", group)) < min_sample_size) next
      if (length(grep("low", group)) < min_sample_size) next
      x <- survival::coxph(survival::Surv(time, event) ~ group, data = expr_clin)
      tmp1 <- broom::tidy(x, exponentiate = T, conf.int = T)
      cox.result[[i]] <- c(gene[i], tmp1)
    }
    res.cox <- data.frame(do.call(rbind, cox.result))
    names(res.cox)[c(1, 3, 5)] <- c("gene", "HR", "Wald_z")

    res_survival[[2]] <- res.cox
    names(res_survival)[[2]] <- "res.cox"

    save(res_survival, file = paste0("output_survival/", project, "_survival_results.rdata"))
    cli::cli_alert_success("Analysis done.")
    return(res_survival)
  }
}
