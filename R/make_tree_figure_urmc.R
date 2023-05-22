# This script is to plot results in standard colors

library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(igraph)
library(ggraph)

# Load trees and metadata
data_dir <- "clean_results/urmc/results_2023-05-15T13:15:38.785201-07:00"
metadata <- read.csv(paste0(data_dir, "/nextstrain_metadata.csv"))

clustertracker_tree <- treeio::read.newick(paste0(data_dir, "/parsnp/parsnp/parsnp.tree"))
per_sample_results <- read.delim(file = paste0(data_dir, "/collated_results_per_sample_with_beast.txt"), sep = " ")

augur_tree <- treeio::read.nhx(paste0(data_dir, "/nextstrain/nextstrain__timetree.nhx"))
augur_node_annotations <- as.data.frame(as_tibble(augur_tree))

beast_tree <- treeio::read.beast(paste0(data_dir, "/beast/combined_chains_mcc.tree"))
beast_node_annotations <- as.data.frame(as_tibble(beast_tree))

# Load cluster data
cluster_data <- read.csv(file = paste0(data_dir, "/hiv_trace/snp_alignment.fasta_user.tn93output.csv"))
distance_threshold <- 0.1

# Define outbreak tips
focal_outbreak_samples <- c(
        "GCF_003952105.1_ASM395210v1_genomic.fna",
        "GCF_003952055.1_ASM395205v1_genomic.fna",
        "GCF_003951125.1_ASM395112v1_genomic.fna",
        "GCF_003950635.1_ASM395063v1_genomic.fna",
        "GCF_003951615.1_ASM395161v1_genomic.fna",
        "GCF_003951535.1_ASM395153v1_genomic.fna",
        "GCF_003950575.1_ASM395057v1_genomic.fna",
        "GCF_003952905.1_ASM395290v1_genomic.fna",
        "GCF_003950595.1_ASM395059v1_genomic.fna",
        "GCF_003951195.1_ASM395119v1_genomic.fna",
        "GCF_003951575.1_ASM395157v1_genomic.fna",
        "GCF_003951635.1_ASM395163v1_genomic.fna",
        "GCF_003951115.1_ASM395111v1_genomic.fna",
        "GCF_003952145.1_ASM395214v1_genomic.fna",
        "GCF_003952035.1_ASM395203v1_genomic.fna")  # These are the CICU outbreak isolates

# Get outbreak clade in each tree
levels_back <- 4

clustertracker_outbreak_mrca <- ape::getMRCA(
    phy = clustertracker_tree, 
    tip = focal_outbreak_samples)
root <- length(clustertracker_tree$tip.label) + 1
node_to_select <- ape::nodepath(phy = clustertracker_tree, from = clustertracker_outbreak_mrca, to = root)[levels_back + 1]
clustertracker_outbreak_tree <- ape::extract.clade(
    phy = clustertracker_tree,
    node = node_to_select)
clustertracker_tipdata <- per_sample_results %>%
    filter(
        method == "clustertracker",
        sample %in% clustertracker_outbreak_tree$tip.label) %>%
    select(sample, introduction_node, region) %>%
    mutate(location_to_plot = case_when(
        sample == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "URMC (unlinked)",
        T ~ region)) %>%
    mutate(location_to_color = case_when(
        sample == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "Focal region",
        sample %in% focal_outbreak_samples ~ "Outbreak",
        T ~ "Other (divisional) region"))

augur_outbreak_mrca <- ape::getMRCA(
    phy = as.phylo(augur_tree),
    tip = focal_outbreak_samples)
augur_outbreak_tree <- treeio::tree_subset(
    augur_tree, 
    node = augur_outbreak_mrca, 
    levels_back = levels_back)
augur_tipdata <- as_tibble(augur_outbreak_tree) %>%
    mutate(location_to_plot = case_when(
        label == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "URMC (unlinked)",
        T ~ location)) %>%
    mutate(location_to_color = case_when(
        label == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "Focal region",
        location == "URMC" ~ "Outbreak",
        T ~ "Other (divisional) region")) %>%
    mutate(location_to_color = factor(
        location_to_color, 
        levels = c("Outbreak", "Focal region", "Other (divisional) region")))

beast_outbreak_mrca <- ape::getMRCA(
    phy = as.phylo(beast_tree), 
    tip = focal_outbreak_samples)
beast_outbreak_tree <- treeio::tree_subset(
    beast_tree, 
    node = beast_outbreak_mrca, 
    levels_back = levels_back)
beast_tipdata <- as_tibble(beast_outbreak_tree) %>%
    mutate(location_to_plot = case_when(
        label == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "URMC (unlinked)",
        T ~ location)) %>%
    mutate(location_to_color = case_when(
        label == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "Focal region",
        location == "URMC" ~ "Outbreak",
        T ~ "Other (global) region"))

# Define standard colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

shared_color_scale <- scale_color_manual(
    breaks = c("Outbreak", "Focal region", "Other (divisional) region", "Other (global) region"),
    values = c("#00BA38", "#619CFF", "#F8766D", "#C77CFF"),
    limits = c("Outbreak", "Focal region", "Other (divisional) region", "Other (global) region")) 
shared_theme <- theme(legend.title = element_blank())
tiplab_size <- 2

# Plot clustertracker tree
clustertracker_plot <- ggtree(tr = clustertracker_outbreak_tree) %<+% clustertracker_tipdata +
    geom_tiplab(aes(label = introduction_node, color = location_to_color), size = tiplab_size) +
    shared_color_scale +
    shared_theme +
    lims(x = c(NA, 0.0025)) +
    geom_treescale(label = "subs/site", fontsize = 2)  # y coord relative to # sequences

# Plot augur tree
augur_plot <- ggtree(
    tr = augur_outbreak_tree, 
    aes(color = location_to_color), 
    mrsd = "2017-06-01", 
    as.Date = T) %<+% augur_tipdata +
    geom_tiplab(aes(label = location_to_plot), size = tiplab_size) +
    shared_color_scale +
    theme_tree2() +
    scale_x_date(limits = c(as.Date("1980-01-01"), as.Date("2032-06-01"))) +
    shared_theme

# Plot beast tree (used arbitrary clock rate of 1)
beast_plot <- ggtree(
    tr = beast_outbreak_tree, 
    aes(color = location_to_color)) %<+% beast_tipdata +
    geom_tiplab(aes(label = location_to_plot), size = tiplab_size) +
    shared_color_scale +
    geom_nodelab(
        aes(label = round(as.numeric(posterior) * 100, digits = 0)), 
        size = tiplab_size,
        hjust = 0) +
    shared_theme +
    lims(x = c(NA, 0.0005)) +
    geom_treescale(label = "subs/site", fontsize = 2)  # y coord relative to # sequences

# Plot HIV-TRACE clusters
cluster_data_to_plot <- cluster_data %>%
    filter(ID1 != ID2) %>%
    mutate(is_outbreak_cluster = ID1 %in% focal_outbreak_samples | ID2 %in% focal_outbreak_samples) %>%
    filter(is_outbreak_cluster) %>%
    select(ID1, ID2)
cluster_matrix <- as.matrix(cluster_data_to_plot)
graph <- igraph::graph_from_edgelist(cluster_matrix, directed = F)
# Color nodes
V(graph)$color <- case_when(
    V(graph)$name %in% focal_outbreak_samples ~ "Outbreak", 
    V(graph)$name == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "Focal region",
    T ~ "Other (divisional) region")
# Label nodes
name_to_label <- metadata$location
names(name_to_label) <- metadata$strain
V(graph)$label <- case_when(
    V(graph)$name %in% focal_outbreak_samples ~ "URMC", 
    V(graph)$name == "GCF_003950615.1_ASM395061v1_genomic.fna" ~ "URMC\n(unlinked)",
    T ~ name_to_label[V(graph)$name])
V(graph)$label.color <- V(graph)$color
V(graph)$label.cex <- 0.8

hivtrace_plot <- ggraph::ggraph(graph, layout = "igraph", algorithm = 'kk') +
    ggraph::geom_edge_link(alpha = 0.2) +
    # ggraph::geom_node_point(aes(color = color), size = 3) +
    ggraph::geom_node_text(aes(color = color, label = label), size = tiplab_size) +
    shared_color_scale +
    theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Plot all trees together
legend_plot <- ggpubr::as_ggplot(ggpubr::get_legend(augur_plot))
plot_list(
    clustertracker_plot + theme(legend.position = "none"), 
    augur_plot + theme(legend.position = "none"), 
    beast_plot + theme(legend.position = "none"),
    hivtrace_plot + theme(legend.position = "none"),
    legend_plot,
    labels = c("A", "B", "C", "D"),
    ncol = 3,
    heights = c(1, 0.5))
ggsave(filename = paste0(data_dir, "/all_trees.png"), width = 6, height = 5.5, units = "in")
