% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getmrnaexpr.R
\name{getmrnaexpr}
\alias{getmrnaexpr}
\title{Get GDC TCGA mRNA/lncRNA expression matrix and clinical information}
\usage{
getmrnaexpr(project)
}
\arguments{
\item{project}{valid TCGA project name(s) from 33 TCGA projects:
\itemize{
\item TCGA-ACC
\item TCGA-BLCA
\item TCGA-BRCA
\item TCGA-CESC
\item TCGA-CHOL
\item TCGA-COAD
\item TCGA-DLBC
\item TCGA-ESCA
\item TCGA-GBM
\item TCGA-HNSC
\item TCGA-KICH
\item TCGA-KIRC
\item TCGA-KIRP
\item TCGA-LAML
\item TCGA-LGG
\item TCGA-LIHC
\item TCGA-LUAD
\item TCGA-LUSC
\item TCGA-MESO
\item TCGA-OV
\item TCGA-PAAD
\item TCGA-PCPG
\item TCGA-PRAD
\item TCGA-READ
\item TCGA-SARC
\item TCGA-SKCM
\item TCGA-STAD
\item TCGA-TGCT
\item TCGA-THCA
\item TCGA-THYM
\item TCGA-UCEC
\item TCGA-UCS
\item TCGA-UVM
}

If you provide more than one TCGA project names, it will combine the data.
This is useful when you try to combine different cancer types, such as
TCGA-COAD and TCGA-READ.}
}
\value{
rdata and csv files which are saved in the directory of "output_mRNA_lncRNA_expr".
\itemize{
\item \strong{TCGA-XXX_SummarizedExperiment.rdata}: SummarizedExperiment object,
all the other files are extracted from this object.
\item \strong{TCGA-XXX_clinicalSE.rdata}: indexed clinical information extracted
from the SummarizedExperiment object.
\item \strong{TCGA-XXX_gene_info.rdata}: gene information, including HGNC Gene
Symbol, Ensembl ID, gene type, eta.
\item \strong{TCGA-XXX_mrna_expr_count.rdata}: mRNA count expression matrix
\item \strong{TCGA-XXX_mrna_expr_tpm.rdata}: mRNA tpm expression matrix
\item \strong{TCGA-XXX_mrna_expr_fpkm.rdata}: mRNA fpkm expression matrix
\item \strong{TCGA-XXX_lncrna_expr_count.rdata}: lncRNA count expression matrix
\item \strong{TCGA-XXX_lncrna_expr_tpm.rdata}: lncRNA tpm expression matrix
\item \strong{TCGA-XXX_lncrna_expr_fpkm.rdata}: lncRNA fpkm expression matrix
}
}
\description{
This function provides a powerful workflow to query, download
and prepare the newest TCGA gene expression quantification data and the
clinical information. All you have to supply is a valid TCGA project name. It
can automatically save six types of expression matrix(mRNAcount/tpm/fpkm,
lncRNA count/tpm/fpkm) and the corresponding clinical information, both in
rdata and csv formats.
}
\details{
In GDC database the clinical data can be retrieved from different sources:
\itemize{
\item indexed clinical: a refined clinical data that is created using the XML files
\item XML files: original source of the data
\item BCR Biotab: tsv files parsed from XML files
}

the clinical data extracted from the SummarizedExperiment object is indexed
clinical, and subtype information from marker TCGA papers are added, these
information have prefix "paper_". more information about clinical data,
please see: \href{https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html}{clinical section} of the \code{TCGAbiolinks} vignette.
}
\references{
Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
analysis of TCGA data.” Nucleic Acids Research.
}
