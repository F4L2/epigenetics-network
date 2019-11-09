rm(list = ls()) #clear environment

load("D:/CODING/epigenetics-network/graphs.RData")
load("D:/CODING/epigenetics-network/main.RData")

load("D:/CODING/epigenetics-network/plot.RData") #if disponible

source("D:/CODING/epigenetics-network/func.R")
library(reshape2)
library(ggplot2)
library(igraph)
library(bnlearn) # Bayesian network learning : contains efficient implementations of several algorithms
library(GENIE3) # GEne Network Inference with Ensemble of trees
library(Rgraphviz)
library(miic) 


#no deviation
plot_score_list(hc_scores)
plot_score_list(RF_scores)
plot_score_list(miic_scores)


#alphas
plot_score_list(pc_scores)
plot_score_list(mmhc_scores)


#nedges GENIE3
rf_genie = RF_GENIE[[1]] #enlever les graph superflu 

  #RF
max_row = nrow(rf_genie)
rf_scores = c()
for( l in seq(1,max_row,by=1)){
  print(l)
  genie3_res = as.matrix(rf_genie[1:l, 1:2]) # columns 1 regulators, 2 targets
  graph = graph_from_edgelist(genie3_res)
  
  pred_mat <- as.data.frame(get.edgelist(graph))
  rf_scores = append(rf_scores, fold_enrichment(pred_mat, ref_mat)[[1]] )
}

#plot with just fold_enr
barplot(rf_scores)



  #ET
et_genie = ET_GENIE
et_scores = c()
for( l in seq(1,max_row,by=1)){
  print(l)
  genie3_res = as.matrix(et_genie[1:l, 1:2]) # columns 1 regulators, 2 targets
  graph = graph_from_edgelist(genie3_res)
  
  pred_mat <- as.data.frame(get.edgelist(graph))
  et_scores = append(et_scores, fold_enrichment(pred_mat, ref_mat)[[1]] )
}
barplot(et_scores)

save.image("D:/CODING/epigenetics-network/plot.RData")



#see optimal scores
rf_opti = which.max( rf_scores )

genie3_res = as.matrix(rf_genie[1:rf_opti, 1:2]) # columns 1 regulators, 2 targets
graph = graph_from_edgelist(genie3_res)

pred_mat <- as.data.frame(get.edgelist(graph))
fold_enrichment(pred_mat, ref_mat)[[1]]




et_opti = which.max( et_scores )
genie3_res = as.matrix(et_genie[1:et_opti, 1:2]) # columns 1 regulators, 2 targets
graph = graph_from_edgelist(genie3_res)

pred_mat <- as.data.frame(get.edgelist(graph))
fold_enrichment(pred_mat, ref_mat)[[1]]




#plot best graphs

#HC
hc_res = HC_GRAPHS[[1]]
graphviz.plot(hc_res)

graph = igraph.from.graphNEL(as.graphNEL(hc_res)) #convert bnlearn graph to igraph
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}


#PC
pc_res = PC_GRAPHS[[17]]
graphviz.plot(pc_res)

graph = igraph.from.graphNEL(as.graphNEL(pc_res)) #convert bnlearn graph to igraph
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}

#MMHC
mmhc_res = MMHC_GRAPHS[[19]]
graphviz.plot(mmhc_res)

graph = igraph.from.graphNEL(as.graphNEL(mmhc_res)) #convert bnlearn graph to igraph
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}

#miic
miic_res = MIIC_GRAPHS[[1]]
miic.plot(miic_res, igraphLayout = layout_nicely)

graph = g= graph_from_adjacency_matrix( miic_res["adjMatrix"][[1]], mode = "directed") #miic to igraph
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}

#GENIE3
genie3_rf_res = as.matrix(rf_genie[1:rf_opti, 1:2])
graph_rf = graph_from_edgelist(genie3_rf_res)
plot(graph_rf)

graph = graph_from_edgelist(genie3_rf_res)
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}


genie3_et_res = as.matrix(et_genie[1:et_opti, 1:2])
graph_et = graph_from_edgelist(genie3_et_res)
plot(graph_et)

graph = graph_from_edgelist(genie3_et_res)
edges = get.edgelist(graph)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}


# again with penalty score

rf_genie = RF_GENIE[[1]] #enlever les graph superflu 

#RF
max_row = nrow(rf_genie)
rf_scores = c()
for( l in seq(1,max_row,by=1)){
  print(l)
  genie3_res = as.matrix(rf_genie[1:l, 1:2]) # columns 1 regulators, 2 targets
  graph = graph_from_edgelist(genie3_res)
  
  pred_mat <- as.data.frame(get.edgelist(graph))
  rf_scores = append(rf_scores, fold_enrichment(pred_mat, ref_mat)[[3]] )
}

#plot with just fold_enr
barplot(rf_scores)



#ET
et_genie = ET_GENIE
et_scores = c()
for( l in seq(1,max_row,by=1)){
  print(l)
  genie3_res = as.matrix(et_genie[1:l, 1:2]) # columns 1 regulators, 2 targets
  graph = graph_from_edgelist(genie3_res)
  
  pred_mat <- as.data.frame(get.edgelist(graph))
  et_scores = append(et_scores, fold_enrichment(pred_mat, ref_mat)[[3]] )
}
barplot(et_scores)




#see optimal scores
rf_opti = which.max( rf_scores )

genie3_rf_res = as.matrix(rf_genie[1:rf_opti, 1:2]) # columns 1 regulators, 2 targets
graph_rf = graph_from_edgelist(genie3_rf_res)
plot(graph_rf)
pred_mat <- as.data.frame(get.edgelist(graph_rf))
fold_enrichment(pred_mat, ref_mat)



et_opti = which.max( et_scores )
genie3_res = as.matrix(et_genie[1:et_opti, 1:2]) # columns 1 regulators, 2 targets
graph_et = graph_from_edgelist(genie3_res)
plot(graph_et)
pred_mat <- as.data.frame(get.edgelist(graph_et))
fold_enrichment(pred_mat, ref_mat)





edges = get.edgelist(graph_rf)
autoregul = 0 
for(i in seq(1, length(edges), by=1)){
  link = edges[i,]
  print(i)
  if( link[1] == link[2] ){
    autoregul = autoregul + 1 
  }
}
