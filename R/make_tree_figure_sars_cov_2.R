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
data_dir <- "clean_results/confa_sensitivity/results_replicate_1"
metadata <- read.delim(paste0(data_dir, "/nextstrain_metadata_with_context.tsv"))

clustertracker_tree <- treeio::read.newick(paste0(data_dir, "/iqtree/tree.nwk"))
per_sample_results <- read.delim(file = paste0(data_dir, "/collated_results_per_sample_with_beast.txt"), sep = " ")

augur_tree <- treeio::read.nhx(paste0(data_dir, "/nextstrain/nextstrain__timetree.nhx"))
augur_node_annotations <- as.data.frame(as_tibble(augur_tree))
augur_tree_with_confidence <- treeio::read.nextstrain.json(paste0(data_dir, "/augur_traits/auspice.json"))

beast_tree <- treeio::read.beast(paste0(data_dir, "/beast/mcc.tree"))
beast_node_annotations <- as.data.frame(as_tibble(beast_tree))

# Load cluster data
cluster_data <- read.csv(file = paste0(data_dir, "/hiv_trace/masked_alignment.fasta_user.tn93output.csv"))
distance_threshold <- 0.1

# Define outbreak tips
focal_outbreak_samples <- c(
        "MT520428.1",
        "MT520429.1",
        "MT520386.1",
        "MT520192.1",
        "MT520178.1",
        "MT520441.1",
        "MT520375.1",
        "MT520325.1",
        "MT520446.1",
        "MT520509.1",
        "MT520408.1",
        "MT520316.1",
        "MT520517.1",
        "MT520503.1",
        "MT520459.1",
        "MT520287.1",
        "MT520191.1",
        "MT520470.1",
        "MT520341.1",
        "MT520327.1",
        "MT520547.1",
        "MT520472.1",
        "MT520396.1",
        "MT520314.1",
        "MT520388.1",
        "MT520320.1",
        "MT520282.1",
        "MT520310.1")

# Get outbreak clade in each tree
clustertracker_outbreak_mrca <- ape::getMRCA(
    phy = clustertracker_tree, 
    tip = focal_outbreak_samples)
clustertracker_outbreak_tree <- ape::extract.clade(
    phy = clustertracker_tree,
    node = clustertracker_outbreak_mrca)
clustertracker_tipdata <- per_sample_results %>%
    filter(
        method == "clustertracker",
        sample %in% clustertracker_outbreak_tree$tip.label) %>%
    select(sample, introduction_node, region) %>%
    mutate(location_to_plot = region) %>%
    mutate(location_to_color = case_when(
        sample %in% focal_outbreak_samples ~ "Outbreak",
        region == "Massachusetts" ~ "Focal region",
        T ~ "Other (divisional) region"))

augur_outbreak_mrca <- ape::getMRCA(
    phy = as.phylo(augur_tree),
    tip = focal_outbreak_samples)
augur_outbreak_tree <- treeio::tree_subset(
    augur_tree, 
    node = augur_outbreak_mrca, 
    levels_back = 0)
augur_tipdata <- as_tibble(augur_outbreak_tree) %>%
    mutate(location_to_plot = division) %>%
    mutate(location_to_color = case_when(
        division == "CONF_A" ~ "Outbreak",
        division == "Massachusetts" ~ "Focal region",
        T ~ "Other (divisional) region"))

# Get confidence for outbreak clade root location
augur_tipdata_all <- as_tibble(augur_tree)
augur_outbreak_mcra_label <- augur_tipdata_all[[augur_outbreak_mrca, "label"]]
augur_nodepie_data <- as_tibble(augur_tree_with_confidence) %>%
    filter(label == augur_outbreak_mcra_label) %>%
    select(node, starts_with("division.confidence.")) %>%
    pivot_longer(cols = starts_with("division.confidence."), names_to = "location", values_to = "confidence") %>%
    mutate(location_to_color = case_when(
        location == "division.confidence.CONF_A" ~ "Outbreak",
        location == "division.confidence.Massachusetts" ~ "Focal region",
        T ~ "Other (divisional) region")) %>%
    group_by(node, location_to_color) %>%
    summarize(confidence = sum(as.numeric(confidence), na.rm = T)) %>%
    pivot_wider(id_cols = node, names_from = location_to_color, values_from = confidence) %>%
    mutate(node = length(as.phylo(augur_outbreak_tree)$tip.label) + 1)  # node numbers don't match up, but node labels should

beast_outbreak_mrca <- ape::getMRCA(
    phy = as.phylo(beast_tree), 
    tip = focal_outbreak_samples)
beast_outbreak_tree <- treeio::tree_subset(
    beast_tree, 
    node = beast_outbreak_mrca, 
    levels_back = 0)
beast_tipdata <- as_tibble(beast_outbreak_tree) %>%
    mutate(location_to_plot = region) %>%
    mutate(location_to_color = case_when(
        region == "CONF_A" ~ "Outbreak",
        T ~ "Other (global) region"))

# A convoluted way to parse the list entries for region.set and region.set.prob
# from beast node data
regions <- unique(as_tibble(beast_tree)$region)
beast_nodepie_data <- beast_tipdata %>%
    select(node, region.set, region.set.prob) %>%
    unnest_wider(col = region.set, names_sep = "_") %>%
    separate(
        col = region.set.prob, 
        into = c(paste0("prob", 1:length(regions))), 
        sep = ",") %>%
    pivot_longer(
        cols = c(paste0("region.set_", 1:length(regions))), 
        names_to = "region_idx", 
        names_prefix = "region.set_",
        values_to = "region_name") %>%
    pivot_longer(
        cols = c(paste0("prob", 1:length(regions))),
        names_to = "prob_idx", 
        names_prefix = "prob",
        values_to = "prob") %>%
    mutate(prob = as.numeric(gsub(x = prob, pattern = "[c()]", replacement = ""))) %>%
    filter(region_idx == prob_idx, !is.na(region_name)) %>%
    mutate(location_to_plot = case_when(
        region_name == "CONF_A" ~ "Outbreak",
        T ~ "Other (global) region")) %>%
    group_by(node, location_to_plot) %>%
    summarize(prob = sum(prob)) %>%
    pivot_wider(id_cols = node, names_from = location_to_plot, values_from = prob)

# Define standard colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

shared_color_scale <- scale_color_manual(
    breaks = c("Outbreak", "Focal region", "Other (divisional) region", "Other (global) region"),
    values = c("#00BA38", "#619CFF", "#F8766D", "#C77CFF"),
    limits = c("Outbreak", "Focal region", "Other (divisional) region", "Other (global) region"),
    aesthetics = c("color", "fill")) 
shared_theme <- theme(legend.title = element_blank())
tiplab_size <- 1.8

# Plot clustertracker tree
clustertracker_plot <- ggtree(tr = clustertracker_outbreak_tree) %<+% clustertracker_tipdata +
    geom_tiplab(aes(label = introduction_node, color = location_to_color), size = tiplab_size) +
    shared_color_scale +
    shared_theme +
    geom_treescale(label = "subs/site", x = 3.3E-4, y = 3, fontsize = 2)  # y coord relative to # sequences

# Plot augur tree
augur_root_pie <- ggtree::nodepie(augur_nodepie_data, cols = 2:(ncol(augur_nodepie_data)))
augur_root_pie <- lapply(
    X = augur_root_pie, 
    FUN = function(g) g + shared_color_scale)

augur_mrsd <- unlist(metadata %>%
    filter(strain %in% as.phylo(augur_outbreak_tree)$tip.label) %>% 
    arrange(desc(date)) %>%
    head(1) %>%
    select(date))
augur_plot <- ggtree(
    tr = augur_outbreak_tree, 
    aes(color = location_to_color), 
    mrsd = as.Date(augur_mrsd), 
    as.Date = T) %<+% augur_tipdata +
    geom_inset(
        insets = augur_root_pie,
        width = 0.15, height = 0.15,
        hjust = 8) +
    geom_tiplab(aes(label = location_to_plot), size = tiplab_size) +
    shared_color_scale +
    theme_tree2() +
    scale_x_date(limits = c(as.Date("2020-01-01"), as.Date("2020-06-01"))) +
    shared_theme

# Plot beast tree
pies <- ggtree::nodepie(beast_nodepie_data, cols = 2:(ncol(beast_nodepie_data)))
n_tips <- length(as.phylo(beast_outbreak_tree)$tip.label)
internal_pies <- pies[(n_tips + 1):length(pies)]
root_pie <- pies[n_tips + 1]
root_pie <- lapply(
    X = root_pie, 
    FUN = function(g) g + shared_color_scale)
beast_mrsd <- unlist(metadata %>%
    filter(strain %in% as.phylo(beast_outbreak_tree)$tip.label) %>% 
    arrange(desc(date)) %>%
    head(1) %>%
    select(date))
beast_plot <- ggtree(
    tr = beast_outbreak_tree, 
    aes(color = location_to_color), 
    mrsd = as.Date(beast_mrsd), 
    as.Date = T) %<+% beast_tipdata +
    geom_inset(
        insets = root_pie,
        width = 0.2, height = 0.2,
        hjust = 8) +
    geom_tiplab(aes(label = location_to_plot), size = tiplab_size) +
    shared_color_scale +
    geom_nodelab(
        aes(label = round(as.numeric(posterior) * 100, digits = 0)), 
        size = tiplab_size,
        hjust = 0) +
    theme_tree2() +
    scale_x_date(limits = c(as.Date("2020-01-01"), as.Date("2020-06-01"))) +
    shared_theme

# Plot HIV-TRACE clusters
cluster_data_to_graph <- cluster_data %>%
    filter(ID1 != ID2) %>%
    select(ID1, ID2)
graph_all <- igraph::graph_from_edgelist(as.matrix(cluster_data_to_graph), directed = F)
clusters <- cluster_walktrap(graph_all)
outbreak_clusters <- data.frame(
    strain = clusters$names,
    group = clusters$membership) %>%
    group_by(group) %>%
    mutate(is_outbreak = any(strain %in% focal_outbreak_samples)) %>%
    ungroup() %>%
    filter(is_outbreak)
cluster_data_to_plot <- cluster_data_to_graph %>%
    filter(ID1 %in% outbreak_clusters$strain)
graph <- igraph::graph_from_edgelist(as.matrix(cluster_data_to_plot), directed = F)
# Color nodes
focal_region_samples <- metadata %>% filter(division == "Massachusetts")
V(graph)$color <- case_when(
    V(graph)$name %in% focal_outbreak_samples ~ "Outbreak", 
    V(graph)$name %in% focal_region_samples$strain ~ "Focal region",
    T ~ "Other (divisional) region")
# Label nodes
name_to_label <- metadata$division
names(name_to_label) <- metadata$strain
V(graph)$label <- case_when(
    V(graph)$name %in% focal_outbreak_samples ~ "Outbreak", 
    T ~ name_to_label[V(graph)$name])
V(graph)$label.color <- "black"
V(graph)$label.cex <- 0.8

hivtrace_plot <- ggraph::ggraph(graph) +
    ggraph::geom_edge_link(alpha = 0.1) +
    ggraph::geom_node_point(aes(color = color), size = 1) +
    # ggraph::geom_node_text(aes(label = label), size = tiplab_size) +
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
    heights = c(1, 0.15))
ggsave(filename = paste0(data_dir, "/all_trees.png"), width = 6, height = 7, units = "in")
