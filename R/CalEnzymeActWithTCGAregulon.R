#' CalEnzymeActWithTCGAregulon
#'
#' @param eset a gene expression matrix with genes (probes) in rows and samples in columns
#' @param TCGAtype Character string indicating cancer type(Def:HCC). Can be one of the following: "ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"
#' @param ... ... Parameters passed to viper
#'
#' @returnA matrix of inferred activity for each regulator gene in the network across all samples
#' @export
#'
#' @examples
CalEnzymeActWithTCGAregulon <- function(eset,
                                        TCGAtype='HCC',...){

  data(list=paste('TCGA_regulon',TCGAtype,sep = '_'))

  EnzymeActivity <- CalEnzymeAct(eset = eset,
                                 regulon=regulon,...)
  EnzymeActivity

}

