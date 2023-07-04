# Network summarization functions for SmCCNet package
# Original script received from Weixuan L
# Script updated by Brenton Graham (last update: 05/21/2023)

require(dplyr)
require(tidyverse)
require(WGCNA)
require(igraph)

# Network summarization function

' Extract Network Summarization Result from Sub-network
#' 
#' Extract summarization scores (the first 3 prinicipal components) for
#' specified network module with given network size. The proteins will be
#' ranked based on PageRank algorithm, then the top k proteins (where k is the
#' specified subnetwork size) will be included into the final subnetwork to
#' generate the summarization score. For the PC score, the correlation with
#' respect to the phenotype of interest will be calculated and stored. In
#' addition, the correlation between individual proteins and phenotype of
#' interest will also be recorded. The final subnetwork adjacency matrix will
#' be stored into the user-specified working directory of interest.
#' 
#' 
#' @param Abar Adjacency matrix of size pxp extracted from the SmCCA step
#' @param CorrMatrix The correlation matrix calculaed based on X1, it should be
#' pxp as well
#' @param SingleOmicsModule Single omics module after hierarchical clustering,
#' should be in list format.
#' @param data the original protein data.
#' @param Pheno the original phenotype data
#' @param folder the working directory where network is intended to be stored
#' @param ModuleIdx the index of the network module that summarization score is
#' intended to be stored
#' @param pheno_name the phenotype name
#' @param mod_size the preferred network size for the given network module,
#' should be an integer from 1 to the largest possible size of the protein
#' network
#' @param damping damping parameter for the pagerank algorithm
#' @param method Either NetSHy'or 'PCA' indicating which summarization method to use
#' @param network_preference Either 'small' or 'large' indicating whether a smaller or larger network is preferred.
#' @param saving_dir Directory where user prefers to store the result
#' @return a file with all local subnetwork information stored in the local designated directory.
#' 
#' @export
#' 
Network_Summarization_PPR_trim <- function(
    Abar, CorrMatrix, data, Pheno, type, ModuleIdx, max_mod_size, damping = 0.9, 
    method = 'NetSHy', network_preference, saving_dir) {
  
  # Trim module by PPR
  net_ppr <- igraph::graph_from_adjacency_matrix(Abar, weighted = TRUE, diag = FALSE, mode = "undirected")
  igraph::set_vertex_attr(net_ppr, "type", index = igraph::V(net_ppr), as.factor(type))
  
  # All parameter setups are based on the recommendation
  ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = damping, options = list(niter = 10^5, eps = 1e-06))
  
  # Obtain ranked names
  rank_names <- names(sort(ranking$vector, decreasing = TRUE))
  rank_value <- sort(ranking$vector, decreasing = TRUE)
  
  # Final subnetwork size
  mod_size <- 0
  
  # Cut mod size to max_mod_size if # nodes exceeds max_mod_size
  if (nrow(Abar) > max_mod_size) {
    mod_size <- max_mod_size
  }
  
  else {mod_size <- nrow(Abar)}
  
  # Obtain pruned adjacency matrix
  newM.node <- which(colnames(Abar) %in% rank_names[1:mod_size])
  sub_type <- type[newM.node]
  M <- Abar[newM.node, newM.node]
  data_topk_subset <- data[, which(colnames(data) %in% rank_names[1:mod_size])]
  
  # Correlation between selected features and the phenotype
  xzy_correlation <- cor(data_topk_subset, Pheno)
  xzy_correlation_test <- rep(0, ncol(data_topk_subset))
  for (i in 1:ncol(data_topk_subset)) {
    xzy_correlation_test[i] <- cor.test(data_topk_subset[,i], Pheno)$p.value
  }
  
  # Store correlatins and p-values in data frame
  omics_correlation_data <- data.frame(
    name = colnames(data_topk_subset), 
    correlation = xzy_correlation, 
    p = xzy_correlation_test)
  
  # Run PCA/NetSHy using pruned data -------------------------------------------
  
  if (method == 'PCA') {
    
    # Run PCA
    pca_summary <- prcomp(data_topk_subset)
    summary_result <- summary(pca_summary)
    
    # Extract the first three PC scores
    pca_score_df <- data.frame(
      pc1 = pca_summary$x[,1], 
      pc2 = pca_summary$x[,2], 
      pc3 = pca_summary$x[,3], 
      y = Pheno)
    
    # PC loading
    pc_loading <- summary_result[["rotation"]]}
  
  else if(method == 'NetSHy') {
    
    # Run NetSHy
    pca_summary <- hybrid_score(data_topk_subset, M, npc = 3, is_alpha = FALSE)
    
    # Extract the first three PC scores
    pca_score_df <- data.frame(
      pc1 = pca_summary[[1]][,1], 
      pc2 = pca_summary[[1]][,2],
      pc3 = pca_summary[[1]][,3], 
      y = Pheno)
    
    # PC loading
    pc_loading <- pca_summary[[3]]} 
  
  # Correlations between PC scores and phenotype
  pc_correlation <- cor(pca_score_df[, 1:3], pca_score_df[,4])
  subnet_corr <- CorrMatrix[newM.node, newM.node]
  pca_importance <- summary_result$importance
  
  # Report subnetwork information
  cat(paste0('The final network size is: ', nrow(M), ' with PC1 correlation w.r.t. phenotype to be: ', round(pc_correlation[1], 3), '\n'))
  
  # Save objects in Rdata file
  save(
    M, pca_score_df, pc_correlation, pc_loading, omics_correlation_data, subnet_corr, 
    pca_importance, file = paste0(saving_dir, "/subnetwork_", ModuleIdx, ".Rdata"))
}


hybrid_score <- function(X, A, is_alpha = TRUE, npc = 1) {
  
  # ----------------------------------------------------------------------------
  # Function Parameters
  # X: data matrix (n, p)
  # A: corresponding adjacency matrix
  # pc_id: PC index
  # ----------------------------------------------------------------------------
  
  g = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Laplacian
  L2 <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE) 
  L2 <- as.matrix(igraph::graph.laplacian(L2, normalized = F))
  
  # TOM
  TOM_A = WGCNA::TOMsimilarity(as.matrix(A), verbose = FALSE)
  
  alpha = igraph::graph.density(g)
  X = scale(X, center = TRUE, scale = TRUE)
  
  # Weighted approach
  if (is_alpha == TRUE) {
    temp  = (1-alpha)*(X %*% L2) + alpha*(X %*% TOM_A)
  }
  
  else { temp = (X %*% L2) }
  
  temp = summary(prcomp(temp))
  h_score = temp$x[,1:npc]
  importance = temp$importance[, 1:npc]
  loading = temp$rotation[, 1:npc]
  return(list(h_score, importance, loading))
  
}
