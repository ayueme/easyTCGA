% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getclinical.R
\name{getclinical}
\alias{getclinical}
\title{Get GDC TCGA clinical data}
\usage{
getclinical(project)
}
\arguments{
\item{project}{valid TCGA project name from 33 TCGA projects:
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
}}
}
\value{
rdata files which are saved in the directory of "output_clinical".
\itemize{
\item \strong{TCGA-XXX_clinicalXML_admin.rdata}: admin data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_drug.rdata}: drug data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_followUp.rdata}: follow up data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_newTumorEvent.rdata}: new tumor event data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_patient.rdata}: general clinical data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_radiation.rdata}: radiation data extracted from XML files.
\item \strong{TCGA-XXX_clinicalXML_stageEvent.rdata}: pathologic stage data extracted from XML files.
\item \strong{TCGA-XXX_clinical_indexed.rdata}: indexed clinical data.
}
}
\description{
Get TCGA clinical data, including indexed clinical data, general information,
survival, tumor stage, drug, radiation, eta.
}
\details{
In GDC database the clinical data can be retrieved from different sources:
\itemize{
\item indexed clinical: a refined clinical data that is created using the XML files
\item XML files: original source of the data
\item BCR Biotab: tsv files parsed from XML files
}

more information about clinical data, please see: \href{https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html}{clinical section}
of the \code{TCGAbiolinks} vignette.

\code{\link[=getclinical]{getclinical()}} can retrieve both the indexed and XML clinical data.
}
\references{
Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D,
Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G,
Noushmehr H (2015). “ TCGAbiolinks: An R/Bioconductor package for integrative
analysis of TCGA data.” Nucleic Acids Research.
}
