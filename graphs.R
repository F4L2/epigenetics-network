rm(list = ls()) #clear environment

library(igraph)
library(Matrix)
library(ggplot2)
library(minet) # Mutual Information NETworks
library(bnlearn) # Bayesian network learning : contains efficient implementations of several algorithms
library(GENIE3) # GEne Network Inference with Ensemble of trees
library(Rgraphviz)
library(miic) # Learning Causal or Non-Causal Graphical Models Using Information Theory


nbCores = detectCores(all.tests = FALSE, logical = TRUE)
data = read.table("D:/CODING/epigenetics-network/data/norm_preProcessed_32nodes.txt", header = T)
TF_names = colnames(data)


#mutual info
#mi_matrix = minet::build.mim(data)



#ref graph
experimental_g <- read_graph("D:/CODING/epigenetics-network/tools/experimental_graph.txt")
experimental_g <- set.vertex.attribute(experimental_g, "name", value= TF_names)
plot(experimental_g)



#hc_res = bnlearn::hc(data)
#graphviz.plot(hc_res)

HC_GRAPHS <- c()
for (i in seq(1, 5, by = 1) ){
  print(i)
  hc_res = bnlearn::hc(data)
  #igraph.from.graphNEL(as.graphNEL(hc_res)) #convert bnlearn graph to igraph
  HC_GRAPHS <- append(HC_GRAPHS, list(hc_res))
}

save.image("D:/CODING/epigenetics-network/5graphs.Rdata")


#constraint based

#pc_res = bnlearn::pc.stable(data, alpha = 0.05)
#graphviz.plot(pc_res)
#pc_graph1 = igraph.from.graphNEL(as.graphNEL(pc_res)) #convert bnlearn graph to igraph

PC_GRAPHS <- c()
for (i in seq(0.00, 1, by = 0.05) ){
  print(i)
  pc_res = bnlearn::pc.stable(data, alpha = i)
  PC_GRAPHS <- append(PC_GRAPHS, list(pc_res))
  
  
  #igraph.from.graphNEL(as.graphNEL(pc_res)) #convert bnlearn graph to igraph
}
save.image("D:/CODING/epigenetics-network/5graphs.Rdata")



    #miic ,constraint based too


#hybrid
#mmhc_res = bnlearn::mmhc(data, restrict.args = list(alpha=0.05))
#graphviz.plot(mmhc_res)

MMHC_GRAPHS <- c()
for (i in seq(0.05, 0.95, by = 0.05) ){
  
  print( c(i,0) )
  mmhc_res = bnlearn::mmhc(data, restrict.args = list(alpha=i))
  MMHC_GRAPHS <- append(MMHC_GRAPHS, list(mmhc_res))
  print( c(i,1) )
  
  #igraph.from.graphNEL(as.graphNEL(mmhc_res)) #convert bnlearn graph to igraph
}
save.image("D:/CODING/epigenetics-network/mmhcgraphs.Rdata")


#feature importance


#RF/all
RF_GENIE <- c()
for (i in seq(1, 5, by = 1) ){
  print(i)
  genie3_weights = GENIE3(t(data), treeMethod= "RF", K="all", verbose=F, nCores = nbCores)
  RF_GENIE <- append(RF_GENIE, list( getLinkList(genie3_weights) ) )
}
save.image("D:/CODING/epigenetics-network/5graphs.Rdata")


#ET/all
ET_GENIE <- c()
for (i in seq(1, 5, by = 1) ){
  print(i)
  genie3_weights = GENIE3(t(data), treeMethod= "ET", K="all", verbose=F, nCores = nbCores)
  ET_GENIE <- append(ET_GENIE, list( getLinkList(genie3_weights) ) )
}

genie3_weights = GENIE3(t(data), treeMethod= "ET", K="all", verbose=F, nCores = nbCores)
ET_GENIE = getLinkList(genie3_weights)
save.image("D:/CODING/epigenetics-network/ETgraph.Rdata")


#save.image("D:/CODING/epigenetics-network/HC_PC_RFGENIE.Rdata")

#####

#TODO: get the right number of rows in linkList [1:nrow(linkList)]

#nRow = floor( sqrt(nrow(data)) )  #heuristic

#genie3_res = as.matrix(linkList[1:nRow, 1:2]) # columns 1 regulators, 2 targets
#genie3_graph = graph_from_edgelist(genie3_res)
#plot(genie3_graph)







