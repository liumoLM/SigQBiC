% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureQBiC.R
\name{SignatureQBiC}
\alias{SignatureQBiC}
\title{Generate GR and LR for a TF-signature pair (Signature-QBiC)}
\usage{
SignatureQBiC(QBiC_score_file_path, pvalue_file_path, sig, plot.path = NULL)
}
\arguments{
\item{QBiC_score_file_path}{All QBiC scores for a universal PBM. The file can be downloaded from
http://qbic.genome.duke.edu/downloads}

\item{pvalue_file_path}{p_values for all twelvemers. The file can be downloaded from
http://qbic.genome.duke.edu/downloads}

\item{sig}{A one-column matrix contains probability of each mutation types
with row.names set to SigQBiC::mut_types.}
}
\value{
A list containing GR and LR.
}
\description{
SignatureQBiC model. This function generates Gain Ratio and Loss Ratio
           for a TF (represented by universal PBM) with
           a mutational signature (or a mutation spectrum)
}
