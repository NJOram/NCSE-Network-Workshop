
## National Centre for Soil Ecology Day - October 30, 2024 ##
## Network Analysis Workshop ##
## Hosted by: Natalie Oram (UvA: n.j.oram@uva.nl) and Doina Bucur (UTwente: d.bucur@utwente.nl) ##

## This script outlines 7 steps to creating co-occurrence networks 
## The networks are built with with Spiec.Easi:  https://github.com/zdk123/SpiecEasi, visualized with igraph: https://r.igraph.org/ and topology determined with various packages (igraph, braingraph). 
## We will use (part of) the published data from van Rijssel et al., 2022: https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16571, namely, the bacterial ASVs in clay soil, from organic or conventionally managed soil. These data have been aggregated to the genus level to reduce the number of nodes (taxa) that are in the network (i.e. we have total number of reads per bacterial genus). In total, we have data from 64 organic farms and 62 conventional farms, and 317 different bacterial genera. 

# Further reading:
# Kurtz et al., 2015. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226
# Chen et al., 2024. https://arxiv.org/pdf/2407.03897
# Chen and Bucur, 2024. https://link.springer.com/chapter/10.1007/978-3-031-57515-0_13 
# Caruso et al., 2022. https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13985
# Gueseva et al., 2022. https://www.sciencedirect.com/science/article/pii/S003807172200061X
# Shan et al., 2023. https://www.nature.com/articles/s41559-023-02021-z 


rm(list=ls())

# Install the packages if necessary. Sometimes SpiecEasi can be a bit moody, but in the newest version of R (4.4.1) it installs fine :-) 

# install_github("zdk123/SpiecEasi")
# install.packages("remotes")
# remotes::install_github("microbiome/microbiome")
# devtools::install_github("briatte/ggnet")

# Libraries
library(tidyverse) 
library(phyloseq)
library(microbiome)

library(devtools)
library(SpiecEasi) # building network
library(igraph) # visualizing network + topology
library(intergraph) # for final network figures
library(network) # for final network figures
library(graph4lg) # for re-building null model ensemble in R
library(ggnet)

# packages for network topology
library(DescTools) # for calculating area under the curve
library(qgraph) # for calculating small worldness
library(brainGraph) # calculating network robustness

library(ggthemes) # for nice plot layouts
library(cowplot) # to easily combine multiple plots


# Step 1a: Load data and make phyloseq object -------------
## Load data
## The data format we need is a phyloseq object. Here, we combine the OTU table, (simplified) taxonomy table, and the sample metadata (all .csv files) into a phyloseq object. If you have the phloseq object, just start from there. 

otu_table<-read.csv("Example data\\ASV.BacteriaClayfiltGenus.rawcounts.csv", header = TRUE, row.names = 1) 
  s<-otu_table%>%
    pivot_longer(cols = everything(), 
                 names_to = "Full_code", 
                 values_to = "value")%>%
    select(Full_code) %>%
    distinct()

sample_data<-read.csv("Example data\\metadata_VRijssel_2022_VitalSoil.csv")%>%
  mutate(Full_code = str_replace_all(Full_code, " ", "."))%>%
  right_join(s, by = "Full_code")%>%
  column_to_rownames(var = "Full_code")

names(sample_data)
levels(as.factor(sample_data$Management))

sample_data%>%
  group_by(Management) %>%
  summarise(count = n())


taxonomy<-read.csv("Example data\\ASV.BacteriaClayfiltGenus.taxonomy.csv", header = TRUE)%>%
  mutate(x = genus)%>%
  column_to_rownames(var = "x")


## Make phyloseq object
otu_table_phyloseq <- otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
sample_data_phyloseq <- sample_data(sample_data)
taxonomy_phyloseq <- tax_table(as.matrix(taxonomy))

ps1 <- phyloseq(otu_table_phyloseq, sample_data_phyloseq, taxonomy_phyloseq)
ps1_organic<-subset_samples(ps1, Management == "Organic")
ps1_conventional<-subset_samples(ps1, Management == "Conventional")

# check
ps1
ps1_organic # should be 64 samples
ps1_conventional # should be 62 samples


# Step 1b: Filter data ------------ ------------

## Not done for this dataset as it is already aggregated to the genus level, leaving 317 taxa (total read count at genus level) (which is fine for the network to handle). This is one option of reducing the amount of nodes in the dataset. However, if you are interested in species-level co-occurrences, this is obviously not the best solution for you. 

### How you filter data is an important question to consider. One option is to remove low-abundance taxa or taxa that only occur in a small proportion of your experimental units. E.g. if taxa are only present in 2% of your 64 bacterial communities from organic farms, it is likely that they don't play an important role in your co-occurrence network (as their chance to 'co-occur' is quite small). How many nodes a network can handle will depend on the dataset (the variation in data). There's no hard and fast number for this, but around 200-800 nodes generally works fine. 



# Step 2: Network construction ----------- ------------
# Below is the code to construct the network. This takes a long time! So we have also provided the constructed network for you. 

# See https://github.com/zdk123/SpiecEasi for more details

# Ingredients
## method: can be 'mb' (Neighborhood selection) or 'glasso' (inverse covariance selection)
# nlambda: tells the model how many different lambda values to try (how many ways to wire each node pair in the network). This essentially determines how complex your network is. A value between 20-120 is fine. Higher nlambda = more computational time and very high nlambda (> 150 or so) risks over fitting. 
# lambda.min.ratio: determines/constrains the range of nlambda values (i.e. 1e-2 means that the smallest lambda value is 1% of the largest)
# pulsar.params: rep.num = number of times the network is re-sampled (bootstrapped), thresh: sets the stability threshold for edges in the network. 0.05 means that an edge must occur in 5% of the bootstrapped samples to be considered in the final model. 
# sel.criterion: model sparseness is inferred using the Stability Approach to Regularization Selection (StARS), which involves random sub-sampling of the data set to find a network with low variability in the selected set of edges

# The goal is to have a stable network >> target stability threshold (getStability) should be very close to 0.05
# to improve the stability score, you can adjust:
#### nlambda
#### lambda.min.ratio
#### sel.criterion


### Organic management -------------
# We will make one network for each management type (organic or conventional) so that we can later compare the network topology. 

# Construct the network (this is computationally heavy/takes some time, so we have also provided the network for you, see 'Load network')
set.seed(210)
network_organic<- spiec.easi(ps1_organic, method='mb',nlambda=120, lambda.min.ratio=1e-2,
                        pulsar.params=list(rep.num=99, thresh=0.05),sel.criterion='bstars')
getStability(network_organic) # Should be close to 0.05

## Save network 
  saveRDS(network_organic, "Example data\\network_organic.rds") 


## Load network 
  network_organic<-readRDS("Example data\\network_organic.rds")
  getStability(network_organic) # 0.04949096: this is acceptable!

# Convert adjacency network to igraph format
ig_network_organic<- adj2igraph(getRefit(network_organic),diag=TRUE,
                                vertex.attr=list(name=taxa_names(ps1_organic)))

  ## Extract edge list and attributes  - we will use these to calculate the null model ensemble
  edgelist_network_organic<-as_edgelist(ig_network_organic, names = FALSE)
  edgelist_network_organic1<-as.data.frame(edgelist_network_organic)
  write.csv(edgelist_network_organic1,"Example data\\Organic network edgelist.csv")
  
  ## Extract adjacency matrix
  adj_network_organic<-as_adjacency_matrix(ig_network_organic, sparse=FALSE)
  adj_network_organic<-as.data.frame(adj_network_organic)
  write.csv(adj_network_organic,"Example data\\Organic network adj matrix.csv")
  

### Conventional management -------------
  
# Construct the network
set.seed(210)
network_conventional<- spiec.easi(ps1_conventional, method='mb',nlambda=100, lambda.min.ratio=1e-2,
                             pulsar.params=list(rep.num=99, thresh=0.05),sel.criterion='bstars')
getStability(network_conventional) # Should be close to 0.05

## Save network
saveRDS(network_conventional, "Example data\\network_conventional.rds")


## Load network 
network_conventional<-readRDS("Example data\\network_conventional.rds")
getStability(network_conventional) #  0.04975339 : this is acceptable!

# Convert adjacency network to igraph format
ig_network_conventional<- adj2igraph(getRefit(network_conventional),diag=TRUE,
                                vertex.attr=list(name=taxa_names(ps1_conventional)))
  
  ## Extract edge list and attributes 
  edgelist_network_conventional<-as_edgelist(ig_network_conventional, names = FALSE)
  edgelist_network_conventional1<-as.data.frame(edgelist_network_conventional)
  write.csv(edgelist_network_conventional1,"Example data\\Conventional network edgelist.csv")
  
  ## Extract adjacency matrix
  adj_network_conventional<-as_adjacency_matrix(ig_network_conventional, sparse=FALSE)
  adj_network_conventional<-as.data.frame(adj_network_conventional)
  write.csv(adj_network_conventional,"Example data\\Conventional network adj matrix.csv")
  
  
# Step 3: Visualize the network ------------- ----------

#### Organic -------------
# Set weights of nodes
## Extract the OTU table from the phyloseq object used in the network 
otu_table <- otu_table(ps1_organic)
otu_df <- as.data.frame(otu_table) # check that it has 64 columns

vsize<- rowMeans(clr(otu_table))+3 # use OTU object here
vsize_df<-as.data.frame(vsize)

am_coord<- igraph::layout_with_fr(ig_network_organic)

# Simple Plots
plot_organic<-plot(ig_network_organic, layout=am_coord, vertex.size = vsize, vertex.label=NA)


#### Conventional -------------
# Set weights of nodes
## Extract the OTU table from the phyloseq object used in the network 
otu_table <- otu_table(ps1_conventional)
otu_df <- as.data.frame(otu_table) # check that it has 62 columns

vsize<- rowMeans(clr(otu_table))+3 # use OTU object here
vsize_df<-as.data.frame(vsize)
am_coord<- igraph::layout_with_fr(ig_network_conventional)

# Simple Plot
plot_conventional<-plot(ig_network_conventional, layout=am_coord, vertex.size=vsize, vertex.label=NA)


# Step 4: Create the Null Model Ensemble (python) and open (back to R) ---------------- ------------
# Create Ensemble -----------
# In python, use the following code to make 999 random networks based on the edge-list of your observed network
        # import numpy as np
        # import networkx as nx
        # import pandas as pd
        # from NEMtropy import UndirectedGraph, matrix_generator
        # from NEMtropy.network_functions import build_adjacency_from_edgelist
        # 
        # df = pd.read_csv('Organic network edgelist.csv') # Change as necessary
        # 
        # G = nx.from_pandas_edgelist(df,source='V1',target='V2')
        # G = G.to_undirected()
        # 
        # # # Transform into adjacency matrix
        # adj_bin = nx.adjacency_matrix(G).toarray()
        # 
        # dseq= np.sum(adj_bin, axis=1)
        # 
        # graph = UndirectedGraph(adj_bin)
        # 
        # graph.solve_tool(model="cm_exp",
        #                  method="newton",
        #                  initial_guess="random")
        # 
        # graph.ensemble_sampler(999, cpu_n=1, output_dir="Organic_null_ensemble/") # Create a new folder for each observed network (e.g. 1 for organic, 1 for conventional)


# Open the ensemble in R -------------
#### Organic ------------------
base_path <- "Example data/"

#### import edge list files
temp = list.files(path = paste0(base_path, "Organic_null_ensemble"), pattern = "*.txt", full.names = TRUE)
read.delim_mod <- function(file_path) {
  read.delim(file_path, header = FALSE)
}

ensemble = lapply(temp, read.delim_mod)

#### read in observed edgelist and adjacency matrix 
edgelist_obs<-read.csv(file = 'Example data\\Organic network edgelist.csv')[,-1]
adj_obs<-read.csv(file = "Example data\\Organic network adj matrix.csv")%>%
  column_to_rownames(var="X")
colnames(adj_obs)<-seq(1:dim(adj_obs)[1])

#### functions
ga.func<-function(mat){
  mat<-graph_from_adjacency_matrix(mat,mode=c("undirected"))
  mat
}

ge.func<-function(mat){
  mat<-graph_from_edgelist(mat,directed=FALSE)
  mat
}

  ## Checks - do the edgelist and adjacency matrix give the same pairs?
  graph.obs<-ga.func(as.matrix(adj_obs)) 
  graph.obs
  graph.obs<-ge.func(as.matrix(edgelist_obs))
  graph.obs


#### Make ensemble matrices edge list, graphs, and adjacency matrix
edge.f<-function(edgelist){
  newlist <- strsplit(edgelist$V1," ")
  newlist <- do.call(rbind, newlist)
  newlist
}

ens.ed.ls<-lapply(ensemble,edge.f)
graph.ens<-lapply(ens.ed.ls,ge.func)

graph.ens[1]
graph.obs<-ge.func(as.matrix(edgelist_obs))


#### Rebuild disconnected nodes into the ensemble adj
#### ensemble adj
get_asj<-function(graph){
  mat<-as.matrix(as_adjacency_matrix(graph))
  mat
}
ens_adj<-lapply(graph.ens,get_asj)
ens_adj[[1]]


#### rebuild
rebuild_f<-function(obs.m,ens.m){
  # obs.m<-get_asj(graph.obs)
  # ens.m<-ens_adj[[2]]
  
  colnames(obs.m)<-seq(1:dim(obs.m)[1])
  
  rm<-matrix(0,dim(obs.m)[1]-dim(ens.m)[1],dim(ens.m)[1])
  ens.m1<-rbind(ens.m,rm)
  # dim(ens.m1)
  cm<-matrix(0,dim(obs.m)[1],dim(obs.m)[1]-dim(ens.m)[1])
  ens.m2<-cbind(ens.m1,cm)
  dim(ens.m2)
  
  
  colvec<-1+as.numeric(colnames(ens.m))
  disc<-setdiff(seq(1:dim(obs.m)[1]),colvec)
  colnames(ens.m2)<-c(colvec,disc)
  
  row.names(ens.m2) <- colnames(ens.m2)
  order <- as.character(sort(as.numeric(colnames(ens.m2))))
  ens.m2 <- reorder_mat(mat = ens.m2, order = order)
  
  ens.m2
  
}

rebuild_f(get_asj(graph.obs),ens_adj[[2]])
ens.list.adj<-lapply(ens_adj,rebuild_f,obs.m=get_asj(graph.obs))

##### Reconvert full ens.adj into igraph object, re-name for quicker loading
ens.list.graph<-lapply(ens.list.adj,ga.func)
ens.list.graph[[1]]
ens.list.graph[[2]]
Organic_graph<-ens.list.graph # this is now your ensemble of null models for the organic management network



#### Conventional ------------------
base_path <- "Example data/"

#### import edge list files
temp = list.files(path = paste0(base_path, "Conventional_null_ensemble"), pattern = "*.txt", full.names = TRUE)
read.delim_mod <- function(file_path) {
  read.delim(file_path, header = FALSE)
}

ensemble = lapply(temp, read.delim_mod)

#### read in observed edgelist and adjacency matrix 
edgelist_obs<-read.csv(file = 'Example data\\Conventional network edgelist.csv')[,-1]
adj_obs<-read.csv(file = "Example data\\Conventional network adj matrix.csv")%>%
  column_to_rownames(var="X")
colnames(adj_obs)<-seq(1:dim(adj_obs)[1])

#### functions
ga.func<-function(mat){
  mat<-graph_from_adjacency_matrix(mat,mode=c("undirected"))
  mat
}

ge.func<-function(mat){
  mat<-graph_from_edgelist(mat,directed=FALSE)
  mat
}

## Checks - do the edgelist and adjacency matrix give the same pairs?
graph.obs<-ga.func(as.matrix(adj_obs)) 
graph.obs
graph.obs<-ge.func(as.matrix(edgelist_obs))
graph.obs


#### Make ensemble matrices edge list, graphs, and adjacency matrix
edge.f<-function(edgelist){
  newlist <- strsplit(edgelist$V1," ")
  newlist <- do.call(rbind, newlist)
  newlist
}

ens.ed.ls<-lapply(ensemble,edge.f)
graph.ens<-lapply(ens.ed.ls,ge.func)

graph.ens[1]
graph.obs<-ge.func(as.matrix(edgelist_obs))


#### Rebuild disconnected nodes into the ensemble adj
#### ensemble adj
get_asj<-function(graph){
  mat<-as.matrix(as_adjacency_matrix(graph))
  mat
}
ens_adj<-lapply(graph.ens,get_asj)
ens_adj[[1]]


#### rebuild
rebuild_f<-function(obs.m,ens.m){
  # obs.m<-get_asj(graph.obs)
  # ens.m<-ens_adj[[2]]
  
  colnames(obs.m)<-seq(1:dim(obs.m)[1])
  
  rm<-matrix(0,dim(obs.m)[1]-dim(ens.m)[1],dim(ens.m)[1])
  ens.m1<-rbind(ens.m,rm)
  # dim(ens.m1)
  cm<-matrix(0,dim(obs.m)[1],dim(obs.m)[1]-dim(ens.m)[1])
  ens.m2<-cbind(ens.m1,cm)
  dim(ens.m2)
  
  
  colvec<-1+as.numeric(colnames(ens.m))
  disc<-setdiff(seq(1:dim(obs.m)[1]),colvec)
  colnames(ens.m2)<-c(colvec,disc)
  
  row.names(ens.m2) <- colnames(ens.m2)
  order <- as.character(sort(as.numeric(colnames(ens.m2))))
  ens.m2 <- reorder_mat(mat = ens.m2, order = order)
  
  ens.m2
  
}

rebuild_f(get_asj(graph.obs),ens_adj[[2]])
ens.list.adj<-lapply(ens_adj,rebuild_f,obs.m=get_asj(graph.obs))

##### Reconvert full ens.adj into igraph object, re-name for quicker loading
ens.list.graph<-lapply(ens.list.adj,ga.func)
ens.list.graph[[1]]
ens.list.graph[[2]]
Conventional_graph<-ens.list.graph # This is now your ensemble of null models for your Conventional network


# Step 5: Calculate network topology ------------------ ----------------

# Calculate the measures of topology that are of interest and bind them into a dataframe 
# There are many, many different measures of topology to choose from. Generally, we are interested in network stability and complexity. 
# Check out igraph for many options. Measures used for neural or social networks may also be interesting (e.g. the robustness function we use here was made to analyze MRI data)


# Organic (observed network)-------------------------- ----------------

#### Measures related to STABILITY -----------------

#### Robustness ---------------------
# Run robustness analysis by removing vertices (nodes) based on degree connectivity
robust <- robustness(ig_network_organic, type = c("vertex"), measure = c("degree"))

# check
ggplot(data=robust, aes(x=removed.pct, y=comp.pct))+
  geom_line(linewidth = 1.5)

#### AUC -------------
#Lmax = maximum number of possible links in network = N(N-1)/2 (definition from network science book)
## N = total number of nodes in the network
network<-ig_network_organic
Lmax <- (vcount(network) * (vcount(network) - 1))/2

auc<-robust%>%
  summarise(AUC = AUC(removed.pct,comp.pct))%>%
  mutate(Lmax = Lmax, 
         AUC_norm = AUC/Lmax)


#### Measures related to COMPLEXITY --------------------

#### Transitivity ------------------------
trans<-transitivity(ig_network_organic, type = "average") 
trans<- trans %>%
  as.data.frame()%>%
  dplyr::rename("transitivity" = ".")
print(trans)


#### Betweenness centrality --------------------
between<-mean(igraph::betweenness(ig_network_organic))
between<-between%>%
  as.data.frame()%>%
  dplyr::rename("betweenness_centrality" = ".")
print(between)


#### Modularity ----------------------
mod = cluster_fast_greedy(ig_network_organic) # use the igraph object where all weights are +
print(mod)
mean_mod<-modularity(mod)
mean_mod<-as.data.frame(mean_mod)
print(mean_mod)


### Module membership ------------------
mod_member<-membership(mod)
mod_member_organic<-as.data.frame(mod_member)%>%
  rename("module" = "x")%>%
  rownames_to_column(var = "ASV")%>%
  mutate(management="organic")


# Combine network properties --------------------
list_data = list(trans, mean_mod, between, auc)
netprop_organic<- list_data %>% 
  purrr::reduce(cbind)%>%
  pivot_longer(cols = c(1:6), names_to = "network_property", values_to = "value")%>%
  mutate(management="organic")
print(netprop_organic)

mod_member_organic<-mod_member

robust_organic<-robust%>%
  mutate(management = "organic")



# Conventional (observed network)------------------------- ---------
#### Measures related to STABILITY --------------------

#### Robustness ---------------------
# Run robustness analysis by removing vertices (nodes) based on degree connectivity
robust <- robustness(ig_network_conventional, type = c("vertex"), measure = c("degree"))

# check
ggplot(data=robust, aes(x=removed.pct, y=comp.pct))+
  geom_line(linewidth = 1.5)

#### AUC -------------
#Lmax = maximum number of possible links in network = N(N-1)/2 (definition from network science book)
## N = total number of nodes in the network
network<-ig_network_conventional
Lmax <- (vcount(network) * (vcount(network) - 1))/2

auc<-robust%>%
  summarise(AUC = AUC(removed.pct,comp.pct))%>%
  mutate(Lmax = Lmax, 
         AUC_norm = AUC/Lmax)



#### Measures related to COMPLEXITY -----------------------

#### Transitivity ------------------------
trans<-transitivity(ig_network_conventional, type = "average") 
trans<- trans %>%
  as.data.frame()%>%
  dplyr::rename("transitivity" = ".")
print(trans)


#### Betweenness centrality --------------------
between<-mean(igraph::betweenness(ig_network_conventional))
between<-between%>%
  as.data.frame()%>%
  dplyr::rename("betweenness_centrality" = ".")
print(between)


#### Modularity ----------------------
mod = cluster_fast_greedy(ig_network_conventional) # use the igraph object where all weights are +
print(mod)
mean_mod<-modularity(mod)
mean_mod<-as.data.frame(mean_mod)
print(mean_mod)


### Module membership ------------------
mod_member<-membership(mod)
mod_member_conventional<-as.data.frame(mod_member)%>%
  rename("module" = "x")%>%
  rownames_to_column(var = "ASV")%>%
  mutate(management="conventional")


# Combine network properties --------------------
list_data = list(trans,mean_mod, between, auc)
netprop_conventional<- list_data %>% 
  purrr::reduce(cbind)%>%
  pivot_longer(cols = c(1:6), names_to = "network_property", values_to = "value")%>%
  mutate(management="conventional")
print(netprop_organic)

mod_member_conventional<-mod_member

robust_conventional<-robust%>%
  mutate(management = "conventional")


# Combine network properties from organic and conventional networks (observed) ------ -----------
net_prop_all<-rbind(netprop_conventional, netprop_organic)
write.csv(net_prop_all, "Example data\\net_prop_obs.csv")

robust_all<-rbind(robust_organic, robust_conventional)
write.csv(robust_all, "Example data\\robust_obs.csv")

mod_mem_all<-rbind(mod_member_conventional, mod_member_organic)
write.csv(mod_mem_all, "Example data\\mod_mem_obs.csv")


# Organic (null model ensemble) ----------------- ---------
# 'Organic_graph' is the list of 999 null models

#### Robustness ----------
robust.func<-function(graph){
  graph.r<-robustness(graph, type = c("vertex"), measure = c("degree"))
  graph.r
}
ens.robust<-lapply(Organic_graph,robust.func)

  #### Save
  df <-  as.data.frame(do.call(rbind, ens.robust))
  df1<-df%>%
    group_by(removed.pct)%>%
    summarise(comp.pct = mean(comp.pct))

robust_all_org <-df%>%
  mutate(management="organic",
         model_type ="null") 

robust_mean_org <-df1%>%
  mutate(management="organic",
         model_type ="null") 

  ### Check 
  ggplot(data=df, aes(x=removed.pct, y=comp.pct))+
    geom_line(linewidth = 1.5, alpha = 0.2)+
    geom_line(data=df1, aes(x = removed.pct, y = comp.pct), colour ="black")

#### AUC ----------------------------------
# ens.robust is the list of 999 null models with robustness calculated
df <- tibble(ID = seq_along(ens.robust)) %>%
  mutate(Objects = map(ens.robust, as_tibble))%>%
  unnest(Objects)
head(df)

auc<-df%>%
  group_by(ID)%>%
  summarise(AUC = AUC(removed.pct,comp.pct))%>%
  rename("value" = "AUC")%>%
  dplyr::select(-"ID")%>%
  mutate(variable = "AUC")
print(auc)

#### Transivity -----------------------
transt.func<-function(graph){
  graph.tr<-transitivity(graph, type = "average")
  graph.tr
}

ens.trans<-lapply(Organic_graph,transt.func)

null.d.metric<-unlist(ens.trans)
trans<-as.data.frame(null.d.metric)%>%
  mutate(variable = "transitivity") %>%
  rename("value" = "null.d.metric")

## Check 
ggplot(data = trans, aes(x =null.d.metric ))+
  geom_histogram(color="black", fill="white", bins = 10)+
  theme_bw()

#### Betweenness centrality --------------------
bet.func<-function(graph){
  graph.b<-mean(igraph::betweenness(graph))
  graph.b
}
ens.bet<-lapply(Organic_graph,bet.func)

#### Save
bet <-  as.data.frame(do.call(rbind, ens.bet))%>%
  mutate(variable = "betweenness_centrality") %>%
  rename("value" = "V1")


#### Modularity ----------------------
mod.func<-function(graph){
  graph.m<-modularity(cluster_fast_greedy(graph))
  graph.m
}
ens.mod<-lapply(Organic_graph,mod.func)

#### Save
mod <-  as.data.frame(do.call(rbind, ens.mod))%>%
  mutate(variable = "modularity") %>%
  rename("value" = "V1")

# Bind null ensemble topology --------------------
net_prop_org_ensemble<-rbind(mod,bet,trans,auc)%>%
  mutate(management = "organic")


# Conventional (null model ensemble) ----------------- ---------
# 'Conventional_graph' is the list of 999 null models

#### Robustness ----------
robust.func<-function(graph){
  graph.r<-robustness(graph, type = c("vertex"), measure = c("degree"))
  graph.r
}
ens.robust<-lapply(Conventional_graph,robust.func)

#### Save
df <-  as.data.frame(do.call(rbind, ens.robust))
df1<-df%>%
  group_by(removed.pct)%>%
  summarise(comp.pct = mean(comp.pct))

robust_all_con <-df%>%
  mutate(management="conventional",
         model_type ="null") 

robust_mean_con <-df1%>%
  mutate(management="conventional",
         model_type ="null") 

### Check 
ggplot(data=df, aes(x=removed.pct, y=comp.pct))+
  geom_line(linewidth = 1.5, alpha = 0.2)+
  geom_line(data=df1, aes(x = removed.pct, y = comp.pct), colour ="black")

#### AUC ----------------------------------
# ens.robust is the list of 999 null models with robustness calculated
df <- tibble(ID = seq_along(ens.robust)) %>%
  mutate(Objects = map(ens.robust, as_tibble))%>%
  unnest(Objects)
head(df)

auc<-df%>%
  group_by(ID)%>%
  summarise(AUC = AUC(removed.pct,comp.pct))%>%
  rename("value" = "AUC")%>%
  dplyr::select(-"ID")%>%
  mutate(variable = "AUC")
print(auc)

#### Transivity -----------------------
transt.func<-function(graph){
  graph.tr<-transitivity(graph, type = "average")
  graph.tr
}

ens.trans<-lapply(Conventional_graph,transt.func)

null.d.metric<-unlist(ens.trans)
trans<-as.data.frame(null.d.metric)%>%
  mutate(variable = "transitivity") %>%
  rename("value" = "null.d.metric")

## Check 
ggplot(data = trans, aes(x =null.d.metric ))+
  geom_histogram(color="black", fill="white", bins = 10)+
  theme_bw()

#### Betweenness centrality --------------------
bet.func<-function(graph){
  graph.b<-mean(igraph::betweenness(graph))
  graph.b
}
ens.bet<-lapply(Conventional_graph,bet.func)

#### Save
bet <-  as.data.frame(do.call(rbind, ens.bet))%>%
  mutate(variable = "betweenness_centrality") %>%
  rename("value" = "V1")


#### Modularity ----------------------
mod.func<-function(graph){
  graph.m<-modularity(cluster_fast_greedy(graph))
  graph.m
}
ens.mod<-lapply(Conventional_graph,mod.func)

#### Save
mod <-  as.data.frame(do.call(rbind, ens.mod))%>%
  mutate(variable = "modularity") %>%
  rename("value" = "V1")

# Bind null ensemble topology (conventional)--------------------
net_prop_con_ensemble<-rbind(mod,bet,trans,auc)%>%
  mutate(management = "conventional")

# Bind all ensemble topology -----------------
ensemble_topology<-rbind(net_prop_org_ensemble,net_prop_con_ensemble)
write.csv(ensemble_topology, "Example data\\ensemble_topology.csv")

robust<-rbind(robust_all_con, robust_all_org)
write.csv(robust, "Example data\\robust_ensemble.csv")

robust_mean<-rbind(robust_mean_org, robust_mean_con)
write.csv(robust_mean, "Example data\\robust_ensemble_mean.csv")

# Step 6: compare network topology of the observed networks and the ensemble ------------------ ---------------

# Network properties related to stability 
robust_obs<-read.csv("Example data\\robust_obs.csv")%>%
  mutate(model_type = "observed")
robust_ens_all<-read.csv("Example data\\robust_ensemble.csv")%>%
  mutate(model_type = "null") # all of the ensemble models
robust_ens_mean<-read.csv("Example data\\robust_ensemble_mean.csv")%>%
  mutate(model_type = "null") # mean of the ensemble models
  
# Network properties related to complexity 
net_prop_obs<-read.csv("Example data\\net_prop_obs.csv")%>%
  mutate(model_type = "observed")
net_prop_ens<-read.csv("Example data\\ensemble_topology.csv")%>%
  rename(network_property = variable)%>%
  mutate(model_type = "null")

head(net_prop_obs)
head(net_prop_ens)

net_prop_all<-rbind(net_prop_obs, net_prop_ens)


#### Robustness figure----------
# Observed models
head(robust_obs)
head(robust_ens_mean)

robust_obs1<-robust_obs%>%
  select(c( "X", "comp.pct","removed.pct","management","model_type"))
     
         
robust_all<-rbind(robust_obs1, robust_ens_mean)
names(robust_all)

##### Robustness
plot <- ggplot(data=robust_all, aes(x=removed.pct, y=comp.pct, colour = management, 
                                    linetype = model_type))+
  geom_line(data=robust_ens_all, aes(x=removed.pct, y=comp.pct, colour = management), alpha=0.1, linewidth=1.5)+ # this is the mean of the ensemble
  geom_line(linewidth = 1) +
  scale_colour_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  scale_linetype_manual(values=c(2,1), name = "Model type")+
  labs(title = "Robustness of bacterial networks in conventional and organically managed soil", 
       y = "Size of largest component",
       x = "Percentage of nodes removed", 
       caption = "Solid lines indicate observed network model\nThe shaded area is the 999 null models")+
  theme(axis.title.y=element_text(size=11, colour = "black", margin=(margin(0,10,0,0))))+
  theme(axis.title.x=element_text(size=11, colour = "black", margin=(margin(0,10,0,0))))+
  theme(axis.text=element_text(size=10, hjust = 1))+
  theme(strip.text.x=element_text(size=12))+
  theme_few() 
plot 

##### Area under the curve 
d<-net_prop_all%>%
  filter(network_property == "AUC" & model_type == "null")
d1<-net_prop_all%>%
  filter(network_property == "AUC" & model_type == "observed")

plot1 <-ggplot()+
  geom_violin(data = d, aes(x = management, y = value, fill = management),alpha = 0.3)+
  geom_point(data = d1, aes(x = management, y = value, colour = management), size=3)+
  scale_colour_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  scale_fill_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  labs(y = "Area",
       title = "Area under the curve")+
  theme_few()+
  theme(strip.text.x=element_text(size=12),
        legend.position = "none")
plot1

fig<-plot_grid(plot, plot1, nrow=1, rel_widths = c(5,3))
fig
ggsave("Network robustness.png", plot = fig, width = 14, height = 7, dpi = 300)

#### Network complexity figure -------------------------
names(net_prop_all)
levels(as.factor(net_prop_all$network_property))

##### Betweenness centrality
d<-net_prop_all%>%
  filter(network_property == "betweenness_centrality" & model_type == "null")
d1<-net_prop_all%>%
  filter(network_property == "betweenness_centrality" & model_type == "observed")

plot1 <-ggplot()+
  geom_violin(data = d, aes(x = management, y = value, fill = management),alpha = 0.3)+
  geom_point(data = d1, aes(x = management, y = value, colour = management), size=3)+
  scale_colour_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  scale_fill_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  labs(y = "Number of shortest paths",
       title = "Betweenness centrality")+
  theme_few()+
  theme(strip.text.x=element_text(size=12))
plot1


##### Modularity
d<-net_prop_all%>%
  filter(network_property == "modularity" & model_type == "null")
d1<-net_prop_all%>%
  filter(network_property == "mean_mod" & model_type == "observed")

plot2 <-ggplot()+
  geom_violin(data = d, aes(x = management, y = value, fill = management),alpha = 0.3)+
  geom_point(data = d1, aes(x = management, y = value, colour = management), size=3)+
  scale_colour_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  scale_fill_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  labs(y = "Modularity (Q)",
       title = "Modularity")+
  theme_few()+
  theme(strip.text.x=element_text(size=12))
plot2


##### Transitivity
d<-net_prop_all%>%
  filter(network_property == "transitivity" & model_type == "null")
d1<-net_prop_all%>%
  filter(network_property == "transitivity" & model_type == "observed")

plot3 <-ggplot()+
  geom_violin(data = d, aes(x = management, y = value, fill = management),alpha = 0.3)+
  geom_point(data = d1, aes(x = management, y = value, colour = management), size=3)+
  scale_colour_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  scale_fill_manual(values=c("steelblue",  "goldenrod"), name="Management")+
  labs(y = "Clustering coefficient",
       title = "Transitivity")+
  theme_few()+
  theme(strip.text.x=element_text(size=12))
plot3


### Combine -----------
legend<-get_legend(plot1)
fig<-plot_grid(plot1+theme(legend.position = "none"),
                          plot2+theme(legend.position = "none"),
                          plot3+theme(legend.position = "none"),
                          nrow = 1)
fig1<-plot_grid(fig, legend, rel_widths=c(1,0.1))
fig1
ggsave("Network properties.png", plot = fig1, width = 12, height = 5, dpi = 300, 
       background = "white")



# Step 7: Prepare the final network figures ---------------- --------------------

# Module information
mod_mem<-read.csv("Example data\\mod_mem_obs.csv")%>%
  pivot_longer(cols = 2:318,names_to = "genus",
               values_to = "module")

mod_mem_organic<-mod_mem%>%
  filter(X == "mod_member_organic")%>%
  mutate(genus_name = genus)%>%
  column_to_rownames(var = "genus")
  
mod_mem_conventional<-mod_mem%>%
  filter(X == "mod_member_conventional")%>%
  mutate(genus_name = genus)%>%
  column_to_rownames(var = "genus")


# Organic network -------------
# Re-load the network made in Spiec.Easi/igraph
network_organic<-readRDS("Example data\\network_organic.rds")
temp <- symBeta(getOptBeta(network_organic), mode="maxabs")
weight <- Matrix::summary(t(temp))[,3]

# re-make network in igraph
ig_network_organic<- adj2igraph(getRefit(network_organic),diag=TRUE,
                                edge.attr=list(weight=weight),
                                vertex.attr=list(name=taxa_names(ps1_organic)))


#Increase visualization
organic_net <- asNetwork(ig_network_organic)
network::set.edge.attribute(organic_net, "color", ifelse(organic_net %e% "weight" > 0, "steelblue", "orange"))


### Add information to determine the largest clusters
t<-mod_mem_organic%>%
  group_by(module)%>%
  count()

mod_mem2<-mod_mem_organic%>%
  inner_join(t, by= "module")%>%
  mutate(module = as.factor(module))%>%
  arrange(desc(n))%>%
  mutate(genus = genus_name)%>%
  column_to_rownames(var = "genus")

# Define the order of levels in the legend
desired_order <- mod_mem2%>%
  select(module)%>%
  distinct()
desired_order<-as.character(desired_order$module)

# Define custom labels for the levels  
t<-mod_mem2%>%
  select(module, n)%>%
  distinct()%>%
  mutate(label = paste("Module", module, "n =", n))%>%
  mutate(alpha = ifelse(n > 19, 0.8, 0.3))

t1<-as.numeric(t$alpha)
t1

t2<-t
t2

# Adjust nodesize
otu_table <- otu_table(ps1_organic)
otu_df <- as.data.frame(otu_table) # check that it has 64 columns
nodesize<- rowMeans(clr(otu_table))+3 

vertex_names <- V(ig_network_organic)$name
names(nodesize) <- vertex_names
organic_net %v% "nodesize" <- nodesize[vertex_names]
 # length(organic_net %v% "nodesize") # check the length


#Add module information (which genus is in which module)
mods<-map_levels(rownames(otu_df), from = "genus_name", to = "module", mod_mem2)
m<-as.data.frame(mods)
organic_net %v% "module" <- mods


#Make network plot
mycolors <- scale_color_manual(values = c("#0072B2",  "#117733","#661100","#6a3d9a","#FFEC8B", "#FFA07A", "#E5F5E0", "grey"),
                               breaks = desired_order,  # The order you want
                               labels = t2$label)

set.seed(123)
plot_organic <- ggnet2(organic_net, 
                       node.color = "module", 
                       alpha = "module",
                       edge.alpha = 0.5,
                       label = F, 
                       node.size = "nodesize", 
                       label.size = 2, 
                       edge.color = "color",
                       mode = "fruchtermanreingold")+ 
  guides(color=guide_legend(title=""), size = "none")+
  mycolors+
  labs(title = "Organic management",
       caption = "total nodes/taxa = 317")+
  theme(plot.title = element_text(hjust = 0.5)) 

plot_organic


# Organic network module membership ---------------
mod_mem_organic<-mod_mem_organic%>%
  rename(genus = genus_name)
mod_mem_df<-as.data.frame(rowSums(otu_df))%>%
  rename(reads = `rowSums(otu_df)`)%>%
  rownames_to_column(var = "genus")%>%
  inner_join(mod_mem_organic, by = "genus")%>%
  filter(reads>500)
names(mod_mem_df)

# Figure
plot_mod_mem_organic<-ggplot(data=mod_mem_df)+
  geom_bar(stat="identity", aes(x = genus,y = reads, fill = genus))+
  facet_wrap(.~module, scales = "free_y")+
  theme_few()+
  labs(title = "Module membership of organically managed soil",
       caption = "Bacterial genera with > 500 reads",
       y = "Read number",
       x = "")+
  theme(plot.title = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 18),
        legend.position = "none")+
  coord_flip()
plot_mod_mem_organic

ggsave("Mod_mem_organic.png", plot_mod_mem_organic, width = 20, height = 12, dpi = 300)

# Conventional networks -------------------------------------

# Module information
mod_mem<-read.csv("Example data\\mod_mem_obs.csv")%>%
  pivot_longer(cols = 2:318,names_to = "genus",
               values_to = "module")

mod_mem_conventional<-mod_mem%>%
  filter(X == "mod_member_conventional")%>%
  mutate(genus1 = genus)%>%
  column_to_rownames(var = "genus")%>%
  rename(genus = genus1)

# Re-load the network made in Spiec.Easi/igraph
network_conventional<-readRDS("Example data\\network_conventional.rds")
temp <- symBeta(getOptBeta(network_conventional), mode="maxabs")
weight <- Matrix::summary(t(temp))[,3]

# re-make network in igraph
ig_network_conventional<- adj2igraph(getRefit(network_conventional),diag=TRUE,
                                edge.attr=list(weight=weight),
                                vertex.attr=list(name=taxa_names(ps1_conventional)))


#Increase visualization
conventional_net <- asNetwork(ig_network_conventional)
network::set.edge.attribute(conventional_net, "color", ifelse(conventional_net %e% "weight" > 0, "steelblue", "orange"))


### Add information to determine the largest clusters
t<-mod_mem_conventional%>%
  group_by(module)%>%
  count()

mod_mem2<-mod_mem_conventional%>%
  inner_join(t, by= "module")%>%
  mutate(module = as.factor(module))%>%
  arrange(desc(n))%>%
  mutate(genus1 = genus)%>%
  column_to_rownames(var = "genus")%>%
  rename(genus = genus1)

# Define the order of levels in the legend
desired_order <- mod_mem2%>%
  select(module)%>%
  distinct()
desired_order<-as.character(desired_order$module)

# Define custom labels for the levels  
t<-mod_mem2%>%
  select(module, n)%>%
  distinct()%>%
  mutate(label = paste("Module", module, "n =", n))%>%
  mutate(alpha = ifelse(n > 19, 0.8, 0.3))

t1<-as.numeric(t$alpha)
t1

t2<-t
t2

# Adjust nodesize
otu_table <- otu_table(ps1_conventional)
otu_df <- as.data.frame(otu_table) # check that it has 64 columns
nodesize<- rowMeans(clr(otu_table))+3 

vertex_names <- V(ig_network_conventional)$name
names(nodesize) <- vertex_names
conventional_net %v% "nodesize" <- nodesize[vertex_names]
# length(conventional_net %v% "nodesize") # check the length


#Add module information (which genus is in which module)
mods<-map_levels(rownames(otu_df), from = "genus", to = "module", mod_mem2)
m<-as.data.frame(mods)
conventional_net %v% "module" <- mods


#Make network plot
mycolors <- scale_color_manual(values = c("#0072B2",  "#117733","#661100","#6a3d9a","#FFEC8B", "#FFA07A", "#E5F5E0", 
                                          "cornflowerblue", "darkorchid", "coral3", "cadetblue"),
                               breaks = desired_order,  # The order you want
                               labels = t2$label)

set.seed(123)
plot_conventional <- ggnet2(conventional_net, 
                       node.color = "module", 
                       alpha = "module",
                       edge.alpha = 0.5,
                       label = F, 
                       node.size = "nodesize", 
                       label.size = 2, 
                       edge.color = "color",
                       mode = "fruchtermanreingold")+ 
  guides(color=guide_legend(title=""), size = "none")+
  mycolors+
  labs(title = "Conventional management",
       caption = "total nodes/taxa = 317")+
  theme(plot.title = element_text(hjust = 0.5)) 

plot_conventional


# Conventional network module membership ---------------
mod_mem_df<-as.data.frame(rowSums(otu_df))%>%
  rename(reads = `rowSums(otu_df)`)%>%
  rownames_to_column(var = "genus")%>%
  inner_join(mod_mem_conventional, by = "genus")%>%
  filter(reads>500)
names(mod_mem_df)

# Figure
plot_mod_mem_conventional<-ggplot(data=mod_mem_df)+
  geom_bar(stat="identity", aes(x = genus,y = reads, fill = genus))+
  facet_wrap(.~module, scales = "free_y")+
  theme_few()+
  labs(title = "Module membership of conventionally managed soil",
       caption = "Bacterial genera with > 500 reads",
       y = "Read number",
       x = "")+
  theme(plot.title = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 18),
        legend.position = "none")+
  coord_flip()
plot_mod_mem_conventional
ggsave("Mod_mem_conventional.png", plot_mod_mem_conventional, width = 20, height = 12, dpi = 300)


# Bind and save figures ------------------
p<-plot_grid(plot_organic, plot_conventional)
p
ggsave("Network plots.png", p, width = 18, height = 8, dpi = 300)




