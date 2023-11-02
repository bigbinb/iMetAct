#' getMetGenes
#' Obtaining metabolism-related genes on the metabolic biological network by restarting the walk random algorithm
#'
#' @param network a two-column data.frame representing metabolism-related biological network. This Networkincludes protein-protein interactions, metabolic enzymes-metabolic substrates, metabolic enzymes-metabolites, metabolites-metabolite receptors, metabolites-transporters
#' @param metabolites a vector of strings representing metabolite molecules in a metabolism-related biological network
#' @param filter.pct a filter threshold (def: 0.2). Filter out genes not related to metabolism.
#' @param gamma restart parameter (def: 0.6)
#' @param tmax maximum number of iterations (steps) (def: 1000)
#' @param eps maximum allowed difference between the computed probabilities at the steady state (def. 1e-10)
#' @param norm whether to normalize the adjacency matrix (def: TRUE)
#'
#' @return a vector
#' @export
#'
#' @examples
#' data(MetNetwork)
#' data(metabolites)
#' getMetGenes(network,
#' metabolites,
#' filter.pct=0.2,
#' gamma = 0.6,
#' tmax = 1000,
#' eps = 1e-10,
#' norm = TRUE)
getMetGenes <- function(network,
                        metabolites,
                        filter.pct=0.2,
                        gamma = 0.6,
                        tmax = 1000,
                        eps = 1e-10,
                        norm = TRUE){
  message("Invert to igraph object...")
  igraph_network <- graph_from_data_frame(network,directed = T)
  ADJ = as_adjacency_matrix(igraph_network)

  message("Calculate the affinity score..")
  AffM=RWR(ADJ%>%as.matrix(),
           metabolites,
           gamma = gamma,
           tmax = tmax,
           eps = eps,
           norm = norm)

  affscore <- AffM$p %>% as.data.frame()
  affscore$genes <- rownames(affscore)

  select_genes <- affscore[affscore$. > quantile(affscore$.,filter.pct),]
  select_genes <- select_genes$genes
  select_genes
}
