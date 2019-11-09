rm(list = ls()) #clear environment
update.packages()


install.packages("igraph")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("bnlearn") # Bayesian network learning : contains efficient implementations of several algorithms


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minet")
BiocManager::install("GENIE3")
BiocManager::install("Rgraphviz")



remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
Sys.setenv(MAKEFLAGS = "-j4") # four cores used
install.packages("rstan", type = "source")

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=native",
    "CXX14 = g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=native",
    file = M, sep = "\n", append = TRUE)



# Install MIIC from remote private repository
if (!require(miic)) install.packages(
  "https://miic.curie.fr/download/miic_mixed.tar.gz", 
  repos = NULL, type = "source"
)


#install.packages("funModeling")
