library(tidyverse)
library(stringdist)
library(igraph)
library(RColorBrewer)
library(viridis)
library(ggraph)

# -------- Paths --------
base_path <- "S:/NETWORKS/rework/Clone Network Input NII FFPE/"
results_path <- "S:/NETWORKS/rework/results/1AA/FFPE/"

# -------- Find input files --------
files <- list.files(base_path, pattern = "\\.network_clones\\.csv$", full.names = TRUE)

# -------- Loop --------
for (input_file in files) {
  
  # Extract "X-Y" from filename
  file_name <- basename(input_file)
  sample_id <- sub("\\.network_clones\\.csv$", "", file_name)  # e.g. "2-1"
  
  message("Processing sample: ", sample_id)
  
  # Define output paths dynamically
  #output_con  <- file.path(results_path, "connected", paste0(sample_id, " con freqs.csv"))
  #output_iso  <- file.path(results_path, "isolated",  paste0(sample_id, " iso freqs.csv"))
  graph_full  <- file.path(results_path, "graphs",    paste0(sample_id, "_graph_full.graphml"))
  graph_conn  <- file.path(results_path, "graphs",    paste0(sample_id, "_graph_connected.graphml"))
  membership_out <- file.path(results_path, "membership", paste0(sample_id, " membership.csv"))
  cluster_info_out <- file.path(results_path, "cluster_info", paste0(sample_id, " clusters_info.csv"))
  clones_data_out  <- file.path(results_path, "clones_data",  paste0(sample_id, " clones_data.csv"))
  graph_values_out <- file.path(results_path, "graph_values", paste0(sample_id, " graph_values.csv"))
  plot_out         <- file.path(results_path, "plots",        paste0(sample_id, " graph.tiff"))
  
  # -------- 1. Read Input --------
  clones_data <- read_csv(input_file, show_col_types = FALSE) %>% column_to_rownames(var = colnames(.)[1])
  cdr3_df <- clones_data %>% select(1)
  cdr3_mat <- as.matrix(cdr3_df)
  cdr3_names <- rownames(cdr3_mat)
  
  # -------- 2. Compute Distance Matrix --------
  distance_matrix <- stringdistmatrix(cdr3_mat, cdr3_mat, method = "lv")
  
  # -------- 3. Build Graph (distance ≤ 2 or ≤1 ) --------
  graph <- make_empty_graph(n = nrow(distance_matrix), directed = FALSE)
  V(graph)$name <- cdr3_names
  
  for (i in 1:(nrow(distance_matrix)-1)) {
    for (j in (i+1):nrow(distance_matrix)) {
      if (distance_matrix[i, j] <= 1) {
        graph <- add_edges(graph, c(i, j))
      }
    }
  }
  
  # -------- 3.2. Scale Node Sizes by Frequency --------
  freq_vec <- clones_data$Frequency
  min_size <- min(freq_vec)
  max_size <- max(freq_vec)
  
  clones_data$norm_size <- (freq_vec - min_size) / (max_size - min_size)
  clones_data$scaled_sizes <- min_size + (max_size - min_size) * (clones_data$norm_size * 250) + 3
  
  # -------- 4. Connected vs Isolated Nodes --------
  degrees <- degree(graph)
  connected <- which(degrees != 0)
  isolated  <- which(degrees == 0)
  
  graph_connected <- delete_vertices(graph, isolated)
  graph_isolated  <- delete_vertices(graph, connected)
  
  # -------- 4.1. Graph Metrics Function --------
  compute_graph_metrics <- function(g) {
    eig <- eigen_centrality(g)$vector
    node_df <- tibble(
      clon_ID = V(g)$name,
      degree = degree(g),
      betweenness = betweenness(g),
      closeness = closeness(g),
      eigenvector = eig
    )
    
    graph_df <- tibble(
      num_vertices = vcount(g),
      num_edges = ecount(g),
      transitivity = transitivity(g),
      graph_diameter = diameter(g),
      edge_density = edge_density(g),
      average_path = ifelse(ecount(g) > 0, mean_distance(g), NA),
      assortativity = assortativity_degree(g)
    )
    
    list(node_df = node_df, graph_df = graph_df)
  }
  
  metrics_all <- compute_graph_metrics(graph)
  metrics_connected <- compute_graph_metrics(graph_connected)
  
  num_isolated <- length(isolated)
  num_total <- vcount(graph)
  num_connected <- num_total - num_isolated
  
  isolated_freq <- signif((num_isolated / num_total) * 100, 5)
  connected_freq <- signif((num_connected / num_total) * 100, 5)
  
  graph_summary <- metrics_all$graph_df %>%
    mutate(
      num_isolated = num_isolated,
      isolated_freq = isolated_freq,
      num_connected = num_connected,
      connected_freq = connected_freq,
      edge_density_connected = metrics_connected$graph_df$edge_density
    )
  
  # Save Graphs
  write_graph(graph, graph_full, format = "graphml")
  write_graph(graph_connected, graph_conn, format = "graphml")
  
  # -------- 5. Frequencies --------
  get_freq_df <- function(graph_sub, clones_data) {
    node_names <- V(graph_sub)$name
    clones_data %>%
      rownames_to_column("ID_clon") %>%
      filter(ID_clon %in% node_names) %>%
      select(ID_clon, Frequency)
  }
  
  con_freqs <- get_freq_df(graph_connected, clones_data)
  iso_freqs <- get_freq_df(graph_isolated, clones_data)
  clones_data$Connected <- as.integer(rownames(clones_data) %in% con_freqs[[1]])
  
  # -------- 7. Cluster Analysis --------
  clusters <- components(graph_connected)
  
  membership_df <- data.frame(
    clon_ID = V(graph_connected)$name,
    membership = clusters$membership,
    stringsAsFactors = FALSE
  ) %>%
    left_join(clones_data %>% rownames_to_column("clon_ID"), by = "clon_ID")
  
  cluster_info <- membership_df %>%
    group_by(membership) %>%
    summarise(
      n_clones = n(),
      fr_cluster = sum(Frequency, na.rm = TRUE)
    )
  
  top_clusters <- cluster_info %>% arrange(desc(n_clones)) %>% slice(1:3)
  top_space    <- cluster_info %>% arrange(desc(fr_cluster)) %>% slice(1:3)
  
  cl_summary <- tibble(
    num_vertices = nrow(cdr3_df),
    num_clusters = clusters$no,
    mean_csize = mean(clusters$csize),
    median_csize = median(clusters$csize),
    csize_1 = top_clusters$n_clones[1],
    cspace_1 = top_clusters$fr_cluster[1],
    csize_2 = top_clusters$n_clones[2],
    cspace_2 = top_clusters$fr_cluster[2],
    csize_3 = top_clusters$n_clones[3],
    cspace_3 = top_clusters$fr_cluster[3]
  )
  
  whole_summary <- merge(graph_summary, cl_summary, by = "num_vertices")
  
  normalized_whole_summary <- whole_summary %>%
    mutate(across(
      where(is.numeric) & !matches(c("num_vertices", "isolated_freq", "connected_freq", "cspace_1", "cspace_2", "cspace_3")),
      ~ .x / num_vertices,
      .names = "norm_{.col}"
    ))
  
  # -------- 9. Save Outputs --------
  write.csv(membership_df, membership_out, row.names = FALSE)
  write.csv(cluster_info, cluster_info_out, row.names = FALSE)
  write.csv(clones_data, clones_data_out, row.names = FALSE)
  write.csv(normalized_whole_summary, graph_values_out, row.names = FALSE)
  
  # -------- 10. Plot --------
  clones_data[clones_data$Connected ==1, ] -> connected_data
  V(graph_connected)$cluster <- as.factor(clusters$membership)
  V(graph_connected)$norm_size <- connected_data$norm_size
  
  n_clusters <- max(clusters$membership) + 1
  
  g <- ggraph(graph_connected, layout = "fr") +
    geom_edge_link(width = 1) +
    geom_node_point(aes(size = norm_size, fill = cluster),
                    shape = 21, stroke = 2) +
    scale_fill_manual(values = viridis::inferno(n_clusters, begin = 0.2, end = 0.9)) +
    scale_size_continuous(range = c(2, 17)) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(filename = plot_out,
         plot = g,
         units = "in",
         width = 1000 / 72,
         height = 1036 / 72,
         dpi = 300)
}



##############################################

library(tidyverse)
library(stringdist)
library(igraph)
library(RColorBrewer)
library(viridis)
library(ggraph)

# -------- Paths --------
base_path <- "S:/NETWORKS/rework/Clone Network Input NII FFPE/"
results_path <- "S:/NETWORKS/rework/results/1AA/FFPE/v2/"

# Ensure output directories exist
dir.create(file.path(results_path, "graphs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_path, "membership"), showWarnings = FALSE)
dir.create(file.path(results_path, "cluster_info"), showWarnings = FALSE)
dir.create(file.path(results_path, "clones_data"), showWarnings = FALSE)
dir.create(file.path(results_path, "graph_values"), showWarnings = FALSE)
dir.create(file.path(results_path, "plots"), showWarnings = FALSE)
dir.create(file.path(results_path, "node_metrics"), showWarnings = FALSE)

# -------- Find input files --------
files <- list.files(base_path, pattern = "\\.network_clones\\.csv$", full.names = TRUE)

# -------- Loop --------
for (input_file in files) {
  
  # Extract "X-Y" from filename
  file_name <- basename(input_file)
  sample_id <- sub("\\.network_clones\\.csv$", "", file_name)
  
  message("Processing sample: ", sample_id)
  
  # Define dynamic outputs
  graph_full         <- file.path(results_path, "graphs", paste0(sample_id, "_graph_full.graphml"))
  graph_conn         <- file.path(results_path, "graphs", paste0(sample_id, "_graph_connected.graphml"))
  membership_out     <- file.path(results_path, "membership", paste0(sample_id, " membership.csv"))
  cluster_info_out   <- file.path(results_path, "cluster_info", paste0(sample_id, " clusters_info.csv"))
  clones_data_out    <- file.path(results_path, "clones_data", paste0(sample_id, " clones_data.csv"))
  graph_values_out   <- file.path(results_path, "graph_values", paste0(sample_id, " graph_values.csv"))
  plot_out           <- file.path(results_path, "plots", paste0(sample_id, " graph.tiff"))
  node_metrics_all_out  <- file.path(results_path, "node_metrics", paste0(sample_id, "_node_metrics_all.csv"))
  node_metrics_conn_out <- file.path(results_path, "node_metrics", paste0(sample_id, "_node_metrics_connected.csv"))
  node_metrics_full_out <- file.path(results_path, "node_metrics", paste0(sample_id, "_node_metrics_full.csv"))
  
  # -------- 1. Read Input --------
  clones_data <- read_csv(input_file, show_col_types = FALSE) %>% column_to_rownames(var = colnames(.)[1])
  cdr3_df <- clones_data %>% select(1)
  cdr3_mat <- as.matrix(cdr3_df)
  cdr3_names <- rownames(cdr3_mat)
  
  # -------- 2. Compute Distance Matrix --------
  distance_matrix <- stringdistmatrix(cdr3_mat, cdr3_mat, method = "lv")
  
  # -------- 3. Build Graph (distance ≤ 2 or ≤1 ) --------
  graph <- make_empty_graph(n = nrow(distance_matrix), directed = FALSE)
  V(graph)$name <- cdr3_names
  
  for (i in 1:(nrow(distance_matrix)-1)) {
    for (j in (i+1):nrow(distance_matrix)) {
      if (distance_matrix[i, j] <= 1) {
        graph <- add_edges(graph, c(i, j))
      }
    }
  }
  
  # -------- 3.2 Scale Node Sizes --------
  freq_vec <- clones_data$Frequency
  min_size <- min(freq_vec)
  max_size <- max(freq_vec)
  
  clones_data$norm_size <- (freq_vec - min_size) / (max_size - min_size)
  clones_data$scaled_sizes <- min_size + (max_size - min_size) * (clones_data$norm_size * 250) + 3
  
  # -------- 4. Connected vs Isolated --------
  degrees <- degree(graph)
  connected <- which(degrees != 0)
  isolated  <- which(degrees == 0)
  
  graph_connected <- delete_vertices(graph, isolated)
  graph_isolated  <- delete_vertices(graph, connected)
  
  # -------- 4.1 Compute Metrics Function --------
  compute_graph_metrics <- function(g) {
    
    deg  <- degree(g)
    btw  <- betweenness(g)
    cls  <- closeness(g)
    eig  <- eigen_centrality(g)$vector
    
    # ----- Clone-count vector for normalization -----
    # (Your input CSV always has Frequency column)
    clone_counts <- clones_data$Frequency[match(V(g)$name, rownames(clones_data))]
    
    # ----- Node-level metrics + NORMALIZED versions -----
    node_df <- tibble(
      clon_ID = V(g)$name,
      degree = deg,
      betweenness = btw,
      closeness = cls,
      eigenvector = eig,
      
      # NORMALIZED PER CLONE
      degree_norm = deg / clone_counts,
      betweenness_norm = btw / clone_counts,
      closeness_norm = cls / clone_counts,
      eigenvector_norm = eig / clone_counts
    )
    
    # ----- Graph-level summaries + NORMALIZED summaries -----
    graph_df <- tibble(
      num_vertices = vcount(g),
      num_edges = ecount(g),
      transitivity = transitivity(g),
      graph_diameter = suppressWarnings(diameter(g)),
      edge_density = edge_density(g),
      average_path = ifelse(ecount(g) > 0, mean_distance(g), NA),
      assortativity = suppressWarnings(assortativity_degree(g)),
      
      # NORMAL (RAW) NODE-LEVEL SUMMARIES
      mean_degree = mean(deg),      median_degree = median(deg),
      mean_betweenness = mean(btw), median_betweenness = median(btw),
      mean_closeness = mean(cls),   median_closeness = median(cls),
      mean_eigenvector = mean(eig), median_eigenvector = median(eig),
      
      # NORMALIZED NODE-LEVEL SUMMARIES (your request)
      mean_degree_norm = mean(deg / clone_counts, na.rm = TRUE),
      median_degree_norm = median(deg / clone_counts, na.rm = TRUE),
      
      mean_betweenness_norm = mean(btw / clone_counts, na.rm = TRUE),
      median_betweenness_norm = median(btw / clone_counts, na.rm = TRUE),
      
      mean_closeness_norm = mean(cls / clone_counts, na.rm = TRUE),
      median_closeness_norm = median(cls / clone_counts, na.rm = TRUE),
      
      mean_eigenvector_norm = mean(eig / clone_counts, na.rm = TRUE),
      median_eigenvector_norm = median(eig / clone_counts, na.rm = TRUE)
    )
    
    list(node_df = node_df, graph_df = graph_df)
  }
  
  
  metrics_all <- compute_graph_metrics(graph)
  metrics_connected <- compute_graph_metrics(graph_connected)
  
  # ---- Save node-level metrics ----
  write.csv(metrics_all$node_df, node_metrics_all_out, row.names = FALSE)
  write.csv(metrics_connected$node_df, node_metrics_conn_out, row.names = FALSE)
  
  # -------- 5. Cluster Analysis --------
  clusters <- components(graph_connected)
  
  membership_df <- data.frame(
    clon_ID = V(graph_connected)$name,
    membership = clusters$membership,
    stringsAsFactors = FALSE
  ) %>%
    left_join(clones_data %>% rownames_to_column("clon_ID"), by = "clon_ID")
  
  # -------- 6. Cluster info --------
  cluster_info <- membership_df %>%
    group_by(membership) %>%
    summarise(
      n_clones = n(),
      fr_cluster = sum(Frequency, na.rm = TRUE)
    )
  
  # -------- 7. Whole Summary with Graph-Level Metrics --------
  num_isolated <- length(isolated)
  num_total <- vcount(graph)
  num_connected <- num_total - num_isolated
  
  graph_summary <- metrics_all$graph_df %>%
    mutate(
      num_isolated = num_isolated,
      isolated_freq = num_isolated / num_total * 100,
      num_connected = num_connected,
      connected_freq = num_connected / num_total * 100,
      connected_edge_density = metrics_connected$graph_df$edge_density
    )
  
  # -------- 8. Save Graphs --------
  write_graph(graph, graph_full, format = "graphml")
  write_graph(graph_connected, graph_conn, format = "graphml")
  
  # -------- 9. Save Clone Data --------
  clones_data$Connected <- as.integer(rownames(clones_data) %in% V(graph_connected)$name)
  write.csv(clones_data, clones_data_out, row.names = FALSE)
  
  # -------- 10. Save cluster info and graph values --------
  write.csv(membership_df, membership_out, row.names = FALSE)
  write.csv(cluster_info, cluster_info_out, row.names = FALSE)
  write.csv(graph_summary, graph_values_out, row.names = FALSE)
  
  # -------- 11. Create unified node metrics file --------
  node_metrics_full <- metrics_all$node_df %>%
    left_join(clones_data %>% rownames_to_column("clon_ID"), by = "clon_ID") %>%
    left_join(membership_df %>% select(clon_ID, membership), by = "clon_ID")
  
  write.csv(node_metrics_full, node_metrics_full_out, row.names = FALSE)
  
  # -------- 12. Plot --------
  connected_data <- clones_data[clones_data$Connected == 1, ]
  
  V(graph_connected)$cluster <- as.factor(clusters$membership)
  V(graph_connected)$norm_size <- connected_data$norm_size
  
  n_clusters <- max(clusters$membership) + 1
  
 # g <- ggraph(graph_connected, layout = "fr") +
  #  geom_edge_link(width = 1) +
  #  geom_node_point(aes(size = norm_size, fill = cluster),
   #                 shape = 21, stroke = 2) +
    #scale_fill_manual(values = viridis::inferno(n_clusters, begin = 0.2, end = 0.9)) +
    #scale_size_continuous(range = c(2, 17)) +
  #  theme_void() +
  #  theme(legend.position = "none")
  
#  ggsave(filename = plot_out,
 #        plot = g,
  #       units = "in",
#         width = 1000 / 72,
#         height = 1036 / 72,
     #    dpi = 300)
}

###############################
second window check


library(tidyverse)
library(stringdist)
library(igraph)
library(RColorBrewer)
library(viridis)
library(ggraph)

# -------- Paths --------
base_path <- "S:/NETWORKS/rework/Clone Network Input NII FFPE/"
results_path <- "S:/NETWORKS/rework/results/1AA/FFPE/"

# -------- Find input files --------
files <- list.files(base_path, pattern = "\\.network_clones\\.csv$", full.names = TRUE)

# -------- Loop --------
for (input_file in files) {
  
  # Extract "X-Y" from filename
  file_name <- basename(input_file)
  sample_id <- sub("\\.network_clones\\.csv$", "", file_name)  # e.g. "2-1"
  
  message("Processing sample: ", sample_id)
  
  # Define output paths dynamically
  #output_con  <- file.path(results_path, "connected", paste0(sample_id, " con freqs.csv"))
  #output_iso  <- file.path(results_path, "isolated",  paste0(sample_id, " iso freqs.csv"))
  graph_full  <- file.path(results_path, "graphs",    paste0(sample_id, "_graph_full.graphml"))
  graph_conn  <- file.path(results_path, "graphs",    paste0(sample_id, "_graph_connected.graphml"))
  membership_out <- file.path(results_path, "membership", paste0(sample_id, " membership.csv"))
  cluster_info_out <- file.path(results_path, "cluster_info", paste0(sample_id, " clusters_info.csv"))
  clones_data_out  <- file.path(results_path, "clones_data",  paste0(sample_id, " clones_data.csv"))
  graph_values_out <- file.path(results_path, "graph_values", paste0(sample_id, " graph_values.csv"))
  plot_out         <- file.path(results_path, "plots",        paste0(sample_id, " graph.tiff"))
  
  # -------- 1. Read Input --------
  clones_data <- read_csv(input_file, show_col_types = FALSE) %>% column_to_rownames(var = colnames(.)[1])
  cdr3_df <- clones_data %>% select(1)
  cdr3_mat <- as.matrix(cdr3_df)
  cdr3_names <- rownames(cdr3_mat)
  
  # -------- 2. Compute Distance Matrix --------
  distance_matrix <- stringdistmatrix(cdr3_mat, cdr3_mat, method = "lv")
  
  # -------- 3. Build Graph (distance ≤ 2 or ≤1 ) --------
  graph <- make_empty_graph(n = nrow(distance_matrix), directed = FALSE)
  V(graph)$name <- cdr3_names
  
  for (i in 1:(nrow(distance_matrix)-1)) {
    for (j in (i+1):nrow(distance_matrix)) {
      if (distance_matrix[i, j] <= 1) {
        graph <- add_edges(graph, c(i, j))
      }
    }
  }
  
  # -------- 3.2. Scale Node Sizes by Frequency --------
  freq_vec <- clones_data$Frequency
  min_size <- min(freq_vec)
  max_size <- max(freq_vec)
  
  clones_data$norm_size <- (freq_vec - min_size) / (max_size - min_size)
  clones_data$scaled_sizes <- min_size + (max_size - min_size) * (clones_data$norm_size * 250) + 3
  
  # -------- 4. Connected vs Isolated Nodes --------
  degrees <- degree(graph)
  connected <- which(degrees != 0)
  isolated  <- which(degrees == 0)
  
  graph_connected <- delete_vertices(graph, isolated)
  graph_isolated  <- delete_vertices(graph, connected)
  
  # -------- 4.1. Graph Metrics Function --------
  compute_graph_metrics <- function(g) {
    eig <- eigen_centrality(g)$vector
    node_df <- tibble(
      clon_ID = V(g)$name,
      degree = degree(g),
      betweenness = betweenness(g),
      closeness = closeness(g),
      eigenvector = eig
    )
    
    graph_df <- tibble(
      num_vertices = vcount(g),
      num_edges = ecount(g),
      transitivity = transitivity(g),
      graph_diameter = diameter(g),
      edge_density = edge_density(g),
      average_path = ifelse(ecount(g) > 0, mean_distance(g), NA),
      assortativity = assortativity_degree(g)
    )
    
    list(node_df = node_df, graph_df = graph_df)
  }
  
  metrics_all <- compute_graph_metrics(graph)
  metrics_connected <- compute_graph_metrics(graph_connected)
  
  num_isolated <- length(isolated)
  num_total <- vcount(graph)
  num_connected <- num_total - num_isolated
  
  isolated_freq <- signif((num_isolated / num_total) * 100, 5)
  connected_freq <- signif((num_connected / num_total) * 100, 5)
  
  graph_summary <- metrics_all$graph_df %>%
    mutate(
      num_isolated = num_isolated,
      isolated_freq = isolated_freq,
      num_connected = num_connected,
      connected_freq = connected_freq,
      edge_density_connected = metrics_connected$graph_df$edge_density
    )
  
  # Save Graphs
  write_graph(graph, graph_full, format = "graphml")
  write_graph(graph_connected, graph_conn, format = "graphml")
  
  # -------- 5. Frequencies --------
  get_freq_df <- function(graph_sub, clones_data) {
    node_names <- V(graph_sub)$name
    clones_data %>%
      rownames_to_column("ID_clon") %>%
      filter(ID_clon %in% node_names) %>%
      select(ID_clon, Frequency)
  }
  
  con_freqs <- get_freq_df(graph_connected, clones_data)
  iso_freqs <- get_freq_df(graph_isolated, clones_data)
  clones_data$Connected <- as.integer(rownames(clones_data) %in% con_freqs[[1]])
  
  # -------- 7. Cluster Analysis --------
  clusters <- components(graph_connected)
  
  membership_df <- data.frame(
    clon_ID = V(graph_connected)$name,
    membership = clusters$membership,
    stringsAsFactors = FALSE
  ) %>%
    left_join(clones_data %>% rownames_to_column("clon_ID"), by = "clon_ID")
  
  cluster_info <- membership_df %>%
    group_by(membership) %>%
    summarise(
      n_clones = n(),
      fr_cluster = sum(Frequency, na.rm = TRUE)
    )
  
  top_clusters <- cluster_info %>% arrange(desc(n_clones)) %>% slice(1:3)
  top_space    <- cluster_info %>% arrange(desc(fr_cluster)) %>% slice(1:3)
  
  cl_summary <- tibble(
    num_vertices = nrow(cdr3_df),
    num_clusters = clusters$no,
    mean_csize = mean(clusters$csize),
    median_csize = median(clusters$csize),
    csize_1 = top_clusters$n_clones[1],
    cspace_1 = top_clusters$fr_cluster[1],
    csize_2 = top_clusters$n_clones[2],
    cspace_2 = top_clusters$fr_cluster[2],
    csize_3 = top_clusters$n_clones[3],
    cspace_3 = top_clusters$fr_cluster[3]
  )
  
  whole_summary <- merge(graph_summary, cl_summary, by = "num_vertices")
  
  normalized_whole_summary <- whole_summary %>%
    mutate(across(
      where(is.numeric) & !matches(c("num_vertices", "isolated_freq", "connected_freq", "cspace_1", "cspace_2", "cspace_3")),
      ~ .x / num_vertices,
      .names = "norm_{.col}"
    ))
  
  # -------- 9. Save Outputs --------
  write.csv(membership_df, membership_out, row.names = FALSE)
  write.csv(cluster_info, cluster_info_out, row.names = FALSE)
  write.csv(clones_data, clones_data_out, row.names = FALSE)
  write.csv(normalized_whole_summary, graph_values_out, row.names = FALSE)
  
  # -------- 10. Plot --------
  clones_data[clones_data$Connected ==1, ] -> connected_data
  V(graph_connected)$cluster <- as.factor(clusters$membership)
  V(graph_connected)$norm_size <- connected_data$norm_size
  
  n_clusters <- max(clusters$membership) + 1
  
  g <- ggraph(graph_connected, layout = "fr") +
    geom_edge_link(width = 1) +
    geom_node_point(aes(size = norm_size, fill = cluster),
                    shape = 21, stroke = 2) +
    scale_fill_manual(values = viridis::inferno(n_clusters, begin = 0.2, end = 0.9)) +
    scale_size_continuous(range = c(2, 17)) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(filename = plot_out,
         plot = g,
         units = "in",
         width = 1000 / 72,
         height = 1036 / 72,
         dpi = 300)
}
