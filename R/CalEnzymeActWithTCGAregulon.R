#' CalEnzymeActWithTCGAregulon
#'
#' Calculate metabolic enzyme activity using pre-computed TCGA regulon networks
#'
#' @param eset a gene expression matrix with genes (probes) in rows and samples in columns
#' @param TCGAtype Character string indicating cancer type(Def:LIHC). Can be one of the following: "ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"
#' @param ... ... Parameters passed to viper
#'
#' @returns matrix of inferred activity for each regulator gene in the network across all samples
#' @export
#' @importFrom utils data
#'
#' @examples
#' data(testExp)
#' enz_activity <- CalEnzymeActWithTCGAregulon(testExp, TCGAtype = "LIHC")
CalEnzymeActWithTCGAregulon <- function(eset,
                                        TCGAtype='LIHC',...){

  data(list=paste("TCGA_regulon","_",TCGAtype,sep = ""))

  EnzymeActivity <- CalEnzymeAct(eset = eset,
                                 regulon=regulon,...)
  EnzymeActivity

}
