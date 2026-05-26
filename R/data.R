#' MetNetwork
#'
#' A metabolism-related biological network integrating protein-protein interactions,
#' metabolic enzyme-substrate, enzyme-metabolite, metabolite-receptor, and
#' metabolite-transporter interactions from KEGG, Reactome, Human-GEM, BRENDA,
#' TRRUST2, and STRING databases.
#'
#' @format A data.frame with two columns representing directed interactions:
#' \describe{
#'   \item{source}{source node}
#'   \item{target}{target node}
#' }
"MetNetwork"

#' metabolites
#'
#' A vector of metabolite molecule names used as seed nodes for the
#' Random Walk with Restart (RWR) algorithm.
#'
#' @format A character vector
"metabolites"

#' MetabolicEnzymes
#'
#' A list of known metabolic enzyme genes.
#'
#' @format A character vector
"MetabolicEnzymes"

#' Metabolism_Re_genes
#'
#' Pre-selected metabolism-related genes identified through network propagation.
#'
#' @format A character vector
"Metabolism_Re_genes"

#' network
#'
#' An ARACNe-inferred transcriptional regulatory network with mutual information
#' scores and p-values. Contains regulator-target interactions with
#' columns: Regulator, Target, MI, pvalue.
#'
#' @format A data.frame with 4 columns:
#' \describe{
#'   \item{Regulator}{regulator gene symbol}
#'   \item{Target}{target gene symbol}
#'   \item{MI}{mutual information score}
#'   \item{pvalue}{p-value of the interaction}
#' }
"network"

#' testExp
#'
#' An example gene expression matrix used for testing and demonstration.
#'
#' @format A numeric matrix with genes in rows and samples in columns
"testExp"

#' regulon
#'
#' A pre-computed regulon object derived from TCGA expression data using
#' ARACNe-AP network inference, used by CalEnzymeActWithTCGAregulon().
#'
#' @format A viper regulon object (list of transcriptional regulatory modules)
"regulon"
