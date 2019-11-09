rm(list = ls()) #clear environment
update.packages()

#source("D:/CODING/epigenetics-network/utils.R")
library(igraph)
library(Matrix)
library(ggplot2)
library(minet) # Mutual Information NETworks
library(bnlearn) # Bayesian network learning : contains efficient implementations of several algorithms
library(GENIE3) # GEne Network Inference with Ensemble of trees
library(Rgraphviz)
library(miic) # Learning Causal or Non-Causal Graphical Models Using Information Theory


nbCores = detectCores(all.tests = FALSE, logical = TRUE)


#reconstruct graph from expression data
data = read.table("D:/CODING/epigenetics-network/data/norm_preProcessed_32nodes.txt", header = T)
TF_names = colnames(data)

#ARACNE

mi_matrix = minet::build.mim(data)

aracne_res = minet::aracne(mi_matrix)
aracne_graph = simplify(graph_from_adjacency_matrix(aracne_res, weighted = T, mode = "directed"))
plot(aracne_graph, layout=layout_nicely)



#score based
hc_res = bnlearn::hc(data)
graphviz.plot(hc_res)

a = hc_res["arcs"][[1]]

g= igraph.from.graphNEL(as.graphNEL(hc_res)) #convert bnlearn graph to igraph
plot(g)




#constraint based
pc_res = bnlearn::pc.stable(data, alpha = 0.05)
graphviz.plot(pc_res)



library(funModeling)
library(dplyr)
d_bins=discretize_get_bins(data=data, n_bins=5)
d_data = discretize_df(data=data, data_bins=d_bins , stringsAsFactors=T)

miic_res = miic(d_data, nThreads=nbCores)   ####discretized
miic.plot(miic_res, igraphLayout = layout_nicely)
dev.copy(png,filename="D:/CODING/epigenetics-network/plots/miic.png")
dev.off ()


g= graph_from_adjacency_matrix( miic_res["adjMatrix"][[1]], mode = "directed") #miic to igraph
plot(g)

e = as_edgelist(g) 
link = e[1,] #row, col
tar = link[2]



#identify miic edges
e = miic_res[["edges"]][1]
links = row.names(e)
links[[1]]

split_links = strsplit(links[[1]], "")[[1]]
for (col in seq(1:length(split_links)) ){
  if( split_links[col] == "1"){
    print( col )
  }
}


#hybrid
mmhc_res = bnlearn::mmhc(data, restrict.args = list(alpha=0.05))
graphviz.plot(mmhc_res)


#feature importance

tMethods = c("RF", "ET") #treemethods
KMethods = c("sqrt", "all")
genie3_weights = GENIE3(t(data), treeMethod= "RF", K="all", verbose=F, nCores = nbCores)
linkList = getLinkList(genie3_weights)

#TODO: get the right number of rows in linkList [1:nrow(linkList)]
nRow = floor( sqrt(nrow(data)) )  #heuristic
#select by value > threshold

genie3_res = as.matrix(linkList[1:nRow, 1:2]) # columns 1 regulators, 2 targets
genie3_graph = graph_from_edgelist(genie3_res)
plot(genie3_graph)


typeof(genie3_graph)
E(genie3_graph)






#####
# Analysis

my_hist = function(data, plot_density=FALSE, transform=NULL, theme=theme_classic, vline=NULL){
  varname = deparse(substitute(data))
  plot = ggplot(data.frame(data=data), aes_string(x=data)) + 
    xlab(varname) +
    geom_histogram(aes(y=..density..), alpha=0.5, fill='lightblue', color='grey') +
    theme()
  if(plot_density) plot = plot + geom_density()
  if(!is.null(transform)) plot = plot + scale_x_continuous(trans=transform)
  if(!is.null(vline)) plot = plot + geom_vline(xintercept = vline, color='orange')
  
  plot
}


g = graph_from_data_frame(data, directed = T)
g = simplify(g, remove.multiple = T)

#X11 ()
plot(g) 
dev.copy(png,filename="D:/CODING/epigenetics-network/plots/fullgraph.png")
dev.off ()


vertices = V(g)
edges = E(g)
adj = as_adj(g) #adjacency matrix
as_adj_list(g) #adjacency matrix as a list

degree(g)[c(1,2,5)] #degree of 1,2 and 5
degrees = table(degree(g))
mean_degree = mean(degrees)
median_degree = median(degrees)

my_hist(degrees) #histogram of degrees
dev.copy(png,filename="/home/alex/Documents/RESYS/projet/plots/degree_histogram.png")
dev.off ()

dd_g = degree_distribution(g)
plot(1:length(dd_g), dd_g, type='h', xlab = 'Degree', ylab = 'Density') #distribution of degrees
dev.copy(png,filename="/home/alex/Documents/RESYS/projet/plots/degree_distrib.png")
dev.off ()

plot(1:length(dd_g), dd_g, type='p', log="xy", xlab = 'Degree', ylab = 'Density') #distribution but log
dev.copy(png,filename="/home/alex/Documents/RESYS/projet/plots/degree_distrib_logscale.png")
dev.off ()

graph_mean_path_length = average.path.length(g)
graph_transity = transitivity(g)


#select biggest component of graph
g_big_comp = decompose(g)[[which.max(components(g)$csize)]]
plot(g_big_comp)
dev.copy(png,filename="/home/alex/Documents/RESYS/projet/plots/biggestcomp_graph.png")
dev.off ()

#use biggest component to generate random graph and compare with them. 

N_nodes <- length(V(g_big_comp))
N_edges <- length(E(g_big_comp))
fit_experimental <- fit_power_law(1+degree(g_big_comp, mode="all"), xmin=3)
max_degree <- max(degree(g_big_comp))






# Random graph with uniform probabilty for every edge
g_random <- erdos.renyi.game(N_nodes, N_edges, type='gnm')

# Scale-free network with edge probabilities proportional to node fitness scores.
# You will have to fit a power law to the observed node distribution first to know
#what parameter (alpha) to use.
g_scalefree <- sample_fitness(N_edges, sample((1:max_degree)^-fit_experimental$alpha, 
                                              N_nodes, replace=TRUE))


par(mfrow=c(1,1))
plot(dd_g, log = 'xy', xlab = 'degree', ylab='p(k)',
     main='Node degree distribution', col='red', type = 'p')
points(degree_distribution(g_random), log = 'xy', col='green', type='p')
points(degree_distribution(g_scalefree), log = 'xy', col='skyblue', type='p')

legend('topright', col=c('red', 'green', 'skyblue'), 
       legend = c('Observed PPI, experimental evidence', 'Random with simple probability',
                  'Random with power law'), lty = 1)

# => closer to random with power law. 












#### Explore further

### Fuzzy cognitive map

#library(FCMapper)

matrix <- as.matrix(data)
concept.names = names(data)
results = nochanges.scenario(matrix,iter=25,concept.names)

graph.fcm(matrix,concept.sizes=results$Equilibrium_value,concept.names)


#red arrow = negative value
#blue arrow = positive value
