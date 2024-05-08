library(stringdist)
library(igraph)
library(tidyverse)
library(RColorBrewer)

#OPEN AA DATA
read.csv("C:/.csv",
         header = TRUE, row.names = 1)  %>% as.matrix() -> aa_matrix

#GRAPH AND METRICS CREATION + PLOTTING
#Distance matrix and graph building
rownames(aa_matrix) -> names_clones
distance_matrix <- stringdistmatrix(aa_matrix, aa_matrix, method = "lv")
as.data.frame(distance_matrix) -> distance_df

graph <- graph.empty(n = nrow(distance_matrix), directed = FALSE)
V(graph)$name <- names_clones
for (i in 1:(nrow(distance_matrix) - 1)) {
  for (j in (i + 1):nrow(distance_matrix)) {
    if (distance_matrix[i, j] <= 2) {
      graph <- add_edges(graph, c(i, j))
    }
  }
}


#OPEN FREQUENCY DATA
read.csv("C:/.csv",
         header = TRUE, row.names = 1) -> frequency_df

min_size <- min(frequency_df$Frequency)
max_size <- max(frequency_df$Frequency)
normalized_sizes <- (frequency_df$Frequency - min(frequency_df$Frequency)) / (max(frequency_df$Frequency) - min(frequency_df$Frequency))
scaled_sizes <- min_size + (max_size - min_size) * (normalized_sizes * 250) + 3

#Set the node sizes based on the scaled frequency values
V(graph)$size <- scaled_sizes 
V(graph)$frecuency <- frequency_df

#group up the vertices
#graph = whole graph
isolated <- which(degree(graph) == 0)
connected <- which(degree(graph) != 0)

#graph2 = connected graph
graph2 <- delete_vertices(graph, isolated)
V(graph2)$name %>% as.data.frame() -> node_names
node_names$. -> row.names(node_names)
#graph3 = isolated graph
graph3 <- delete_vertices(graph, connected)
V(graph3)$name %>% as.data.frame() -> iso_names
iso_names$. -> row.names(iso_names)

merge(node_names, frequency_df, by = 0, all = FALSE) -> con_freqs
row.names(con_freqs) <- con_freqs$Row.names
select(con_freqs, Frequency) -> con_freqs

merge(iso_names, frequency_df, by = 0, all = FALSE) -> iso_freqs
row.names(iso_freqs) <- iso_freqs$Row.names
select(iso_freqs, Frequency) -> iso_freqs

as.numeric(length(isolated)) -> num_isolated
as.numeric(length(distance_df)) -> num_total_clones
num_total_clones - num_isolated -> num_connected

(num_isolated/num_total_clones)*100 -> isolated_freq 
signif(isolated_freq, digits = 5) -> isolated_freq
(num_connected/num_total_clones)*100 -> connected_freq 
signif(connected_freq, digits = 5) -> connected_freq

frequency_df[rownames(node_names), colnames(frequency_df)] %>% as.data.frame() -> con_freq
row.names(con_freq) <- row.names(node_names)
colnames(con_freq) <- "frequency"


#clusters
clusters <- components(graph2)

lapply(seq_along(clusters$csize)[clusters$csize>1], function(x)
  V(graph2)$name[clusters$membership %in% x]) -> members
plyr::ldply(members, rbind) %>% t() -> mat_membership

cluster_size <- data.frame(
  row.names = seq.int(length(clusters$csize)),
  n_clones = clusters$csize
)

cl_df <- data.frame(
  row.names =  row.names(node_names),
  clon_ID = row.names(node_names),
  membership = clusters$membership,
  fr_clon = con_freq$frequency
)

cl_df %>% group_by(membership) %>% summarise(fr_cluster = sum(fr_clon)) -> freq_cluster

cbind(cluster_size,
      freq_cluster)-> clusters_info

#extract frequency of a certain cluster 
cl_df[cl_df$membership == 36,] -> nan 
sum(nan$fr_clon)

#coloring by cluster
colorRampPalette(brewer.pal(11, "Spectral"))(50) -> expanded_esp 
V(graph2)$color <- clusters$membership
colors <- brewer.pal(max(clusters$membership) + 1, "Spectral")[clusters$membership + 1]
expanded_esp[clusters$membership + 1] -> colors
V(graph2)$name <- ""
plot(graph2)

#GRAPH METRICS
#node and edge
evcent(graph) -> eigenvector2
node_df <- data.frame(
  row.names = names_clones,
  degree = degree(graph),
  betweenness = betweenness(graph),
  closeness = closeness(graph), 
  ev_centrality = eigenvector2$vector
)

graph_df2 <- data.frame(
  num_vertices = vcount(graph),
  num_edges = ecount(graph),
  transitivity = transitivity(graph),
  graph_diameter = diameter(graph),
  graph_density = graph.density(graph),
  average_path = average.path.length(graph),
  assortativity = assortativity_degree(graph),
  num_isolated = num_isolated,
  isolated_freq = isolated_freq,
  num_connected = num_connected,
  connected_freq = connected_freq,
  graph_density_con = graph.density(graph2)
)

evcent(graph2) -> eigenvector3
node_df2 <- data.frame(
  row.names = node_names$.,
  degree = degree(graph2),
  betweenness = betweenness(graph2),
  closeness = closeness(graph2), 
  ev_centrality = eigenvector3$vector
)

#CLUSTER METRICS
mean(clusters$csize) -> mean_csize
median(clusters$csize) -> median_csize

#top 3 clo size 
top3 <- clusters_info %>%
  arrange(desc(n_clones)) %>%
  #group_by(membership) %>%
  slice(1:3)

top3space <- clusters_info %>%
  arrange(desc(fr_cluster)) %>%
  #group_by(membership) %>%
  slice(1:3)


data.frame(
  num_cl = clusters$no,
  mean_csize = mean_csize,
  median_csize = median_csize,
  csize_1 = top3$n_clones[1],
  csize_2 = top3$n_clones[2],
  csize_3 = top3$n_clones[3],
  cspace_1 = top3$fr_cluster[1],
  cspace_2 = top3$fr_cluster[2],
  cspace_3 = top3$fr_cluster[3]
) -> clnumber

data.frame(
  csize_1_sp = top3space$n_clones[1],
  csize_2_sp = top3space$n_clones[2],
  csize_3_sp = top3space$n_clones[3],
  cspace_1_sp = top3space$fr_cluster[1],
  cspace_2_sp = top3space$fr_cluster[2],
  cspace_3_sp = top3space$fr_cluster[3]
) -> clspace
merge(clnumber, clspace)-> pc_clusters

#save all metrics' dataframes, structured by information level
write.csv(node_df, "/.csv", row.names = TRUE, col.names = TRUE)
write.csv(graph_df2, "/.csv", row.names = TRUE, col.names = TRUE)
write.csv(node_df2, "/.csv", row.names = TRUE, col.names = TRUE)

write.csv(con_freqs, "/.csv", row.names = TRUE, col.names = TRUE)
write.csv(iso_freqs, "/.csv", row.names = TRUE, col.names = TRUE)

write.csv(mat_membership, "/.csv", row.names = TRUE, col.names = TRUE)
write.csv(clusters_info, "/.csv", row.names = TRUE, col.names = TRUE)
write.csv(pc_clusters, "/.csv", row.names = TRUE, col.names = TRUE)
