rm(list = ls()) #clear environment

library(igraph)
library(Matrix)
library(ggplot2)
library(minet) # Mutual Information NETworks
library(bnlearn) # Bayesian network learning : contains efficient implementations of several algorithms
library(GENIE3) # GEne Network Inference with Ensemble of trees
library(Rgraphviz)
library(miic) # Learning Causal or Non-Causal Graphical Models Using Information Theory

source("D:/CODING/epigenetics-network/func.R")
load("D:/CODING/epigenetics-network/graphs.RData")



#ref graph
experimental_g <- read_graph("D:/CODING/epigenetics-network/tools/experimental_graph.txt")
experimental_g <- set.vertex.attribute(experimental_g, "name", value= TF_names)
plot(experimental_g)
ref_mat <- as.data.frame(get.edgelist(experimental_g))


hc_scores = c()
for(i in seq(1, length(HC_GRAPHS), by = 1) ){
  g = HC_GRAPHS[[i]]
  graph = igraph.from.graphNEL(as.graphNEL(g)) #convert bnlearn graph to igraph
  pred_mat <- as.data.frame(get.edgelist(graph))
  hc_scores = append(hc_scores, list(fold_enrichment(pred_mat, ref_mat)) )
}

pc_scores = c()
for(i in seq(2, length(PC_GRAPHS), by = 1) ){
  g = PC_GRAPHS[[i]]
  graph = igraph.from.graphNEL(as.graphNEL(g)) #convert bnlearn graph to igraph
  pred_mat <- as.data.frame(get.edgelist(graph))
  pc_scores = append(pc_scores, list(fold_enrichment(pred_mat, ref_mat)) )
}

mmhc_scores = c()
for(i in seq(2, length(MMHC_GRAPHS), by = 1) ){
  g = MMHC_GRAPHS[[i]]
  graph = igraph.from.graphNEL(as.graphNEL(g)) #convert bnlearn graph to igraph
  pred_mat <- as.data.frame(get.edgelist(graph))
  mmhc_scores = append(mmhc_scores, list(fold_enrichment(pred_mat, ref_mat)) )
}


miic_scores = c()
for(i in seq(1, length(MIIC_GRAPHS), by = 1) ){
  miic_res = MIIC_GRAPHS[[i]]
  graph = g= graph_from_adjacency_matrix( miic_res["adjMatrix"][[1]], mode = "directed") #miic to igraph
  pred_mat <- as.data.frame(get.edgelist(graph))
  miic_scores = append(miic_scores, list(fold_enrichment(pred_mat, ref_mat)) )
}



nRow = floor( sqrt(nrow(data)) )  #heuristic
RF_scores = c()
for(i in seq(1, length(RF_GENIE), by = 1) ){
  linkList = RF_GENIE[[i]]
  genie3_res = as.matrix(linkList[1:nRow, 1:2]) # columns 1 regulators, 2 targets
  graph = graph_from_edgelist(genie3_res)
  
  pred_mat <- as.data.frame(get.edgelist(graph))
  RF_scores = append(RF_scores, list(fold_enrichment(pred_mat, ref_mat)) )
}



linkList = ET_GENIE
genie3_res = as.matrix(linkList[1:nRow, 1:2]) # columns 1 regulators, 2 targets
graph = graph_from_edgelist(genie3_res)

pred_mat <- as.data.frame(get.edgelist(graph))
ET_score = fold_enrichment(pred_mat, ref_mat)




save.image("D:/CODING/epigenetics-network/main.Rdata")
