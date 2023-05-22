# This script is to collate output files from different methods into a unified results file for each run

library(tidyverse)
library(treeio)
library(rjson)
library(ape)

# Analysis-specific parameters
run_names <- c(
    "main",
    "new_tree",
    "replicate_1",
    "replicate_2",
    "replicate_3")
outbreak_annotation <- "CONF_A"

for (run_name in run_names) {
    data_dir <- paste0("clean_results/confa_sensitivity/results_", run_name)
    full_metadata <- read.delim(paste0(data_dir, "/metadata_with_context.tsv"))
    outbreak_metadata <- full_metadata %>% filter(division == outbreak_annotation)
    location_colname <- "division"

    ### Load tabular results files ###

    clustertracker_results <- read.delim(paste0(data_dir, "/clustertracker/introductions.tsv"))


    ### Get clusters from treetime mugration output ###

    mugration_treeio <- treeio::read.nhx(paste0(data_dir, "/nextstrain/nextstrain__timetree.nhx"))
    node_annotations <- as.data.frame(as_tibble(mugration_treeio))
    mugration_phylo <- ape::as.phylo(mugration_treeio)

    # pre-order tree traversal to find outbreak MRCA nodes
    mugration_phylo <- ape::reorder.phylo(mugration_phylo, order = "postorder", index.only = F)
    postorder_edges <- mugration_phylo$edge

    # Traverse tree, finding import events
    is_first <- T
    for (i in nrow(postorder_edges):1) {
        parent <- postorder_edges[i, 1]
        child <- postorder_edges[i, 2]
        parent_annotation <- node_annotations[parent, location_colname]
        child_annotation <- node_annotations[child, location_colname]
        # print(paste("parent:", parent_annotation, "& child:", child_annotation))
        if (is.na(parent_annotation) || child_annotation != parent_annotation) {
            augur_introduction_data_tmp <- data.frame(
                parent = parent,
                child = child,
                parent_annotation = parent_annotation,
                child_annotation = child_annotation
            )
            if (is_first) {
                augur_introduction_data <- augur_introduction_data_tmp
                is_first <- F
            } else {
                augur_introduction_data <- rbind(augur_introduction_data, augur_introduction_data_tmp)
            }
        }  
    }


    ### Get mugration imports in a per-sample tabular format ###

    # function to see if a node is the child of an introduction
    get_introduction <- function(node_idx, augur_introduction_data) {
        if (node_idx %in% augur_introduction_data$child) {
            return(node_idx)
        } else {
            parent <- postorder_edges[which(postorder_edges[, 2] == node_idx), 1]
            return(get_introduction(parent, augur_introduction_data))
        }
    }

    n_tips <- length(mugration_phylo$tip.label)
    is_first <- T
    for (i in 1:n_tips) {
        sample_data_tmp <- data.frame(
            sample = mugration_phylo$tip.label[i],
            region = node_annotations[i, location_colname],
            introduction_node = get_introduction(i, augur_introduction_data))
        if (is_first) {
            sample_data <- sample_data_tmp
            is_first <- F
        } else {
            sample_data <- rbind(sample_data, sample_data_tmp)
        }
    }

    mugration_results <- sample_data %>%
        group_by(introduction_node) %>%
        mutate(cluster_size = n())


    ### Parse HIV-TRACE json output ###

    clusters_json <- readr::read_file(paste0(data_dir, "/hiv_trace/clusters_skip.json"))
    clusters_list <- rjson::fromJSON(json_str = clusters_json)

    cluster_idxs <- clusters_list$trace_results$Nodes$cluster$values
    clustered_sample_names <- clusters_list$trace_results$Nodes$id

    cluster_results <- data.frame(
        sample = clustered_sample_names,
        introduction_node = cluster_idxs)

    outbreak_sample_names <- outbreak_metadata$strain
    sample_names <- full_metadata$strain
    n_clusters <- max(cluster_idxs) + 1  # hivtrace clusters are zero-indexed
    singleton_sample_names <- sample_names[!(sample_names %in% clustered_sample_names)]
    n_singletons <- length(singleton_sample_names)

    singleton_results <- data.frame(
        sample = singleton_sample_names,
        introduction_node = n_clusters:(n_clusters + n_singletons - 1))  # hivtrace clusters are zero-indexed

    hivtrace_results <- rbind(cluster_results, singleton_results) %>% 
        group_by(introduction_node) %>% 
        mutate(cluster_size = n()) %>%
        mutate(region = case_when(
            sample %in% outbreak_sample_names ~ outbreak_annotation, T ~ "Other"))


    ### Combine all per-sample results ###

    full_results <- rbind(
        clustertracker_results %>% 
            dplyr::select(sample, introduction_node, cluster_size, region) %>%
            mutate(method = "clustertracker"),
        mugration_results %>% mutate(method = "mugration"),
        hivtrace_results %>% mutate(method = "hivtrace")
    )

    # inspect results
    tmp <- full_results %>%
        group_by(method, introduction_node) %>%
        summarize(
            n_outbreak_samples = sum(region == outbreak_annotation),
            n_total_samples = n()) %>%
        filter(n_outbreak_samples > 0) %>%
        arrange(method, desc(n_outbreak_samples))
    # print(tmp)

    # save results
    write.table(
        x = full_results, 
        file = paste0(data_dir, "/collated_results_per_sample.txt"),
        row.names = F)

    ### Get per-introduction results ###

    # function to see if a node is the child of an outbreak introduction
    get_outbreak_introduction <- function(node_idx, root_idx, augur_outbreak_introductions) {
        if (node_idx %in% augur_outbreak_introductions$child) {
            # base case - reached outbreak introduction
            return(node_idx)
        } else if (node_idx == root_idx) {
            # base case - reached root and no outbreak introduction found
            return(NA)
        } else {
            parent <- postorder_edges[which(postorder_edges[, 2] == node_idx), 1]
            return(get_outbreak_introduction(parent, root_idx, augur_outbreak_introductions))
        }
    }

    augur_outbreak_introductions <- augur_introduction_data %>% filter(child_annotation == outbreak_annotation)
    n_tips <- length(mugration_phylo$tip.label)
    is_first <- T
    for (i in 1:n_tips) {
        mugration_descendent_data_tmp <- data.frame(
            tip = i,
            label = mugration_phylo$tip.label[i],
            location = node_annotations[i, location_colname],
            introduction_node = get_outbreak_introduction(i, n_tips + 1, augur_outbreak_introductions))
        if (is_first) {
            mugration_descendent_data <- mugration_descendent_data_tmp
            is_first <- F
        } else {
            mugration_descendent_data <- rbind(mugration_descendent_data, mugration_descendent_data_tmp)
        }
    }
    mugration_descendent_data <- mugration_descendent_data %>% 
        filter(!is.na(introduction_node)) %>%
        mutate(method = "mugration")

    ### Get clusters from clustertracker output ###

    full_metadata_tmp <- full_metadata  %>%
        rename("location" = all_of(location_colname))
    node_to_leaves <- read.csv(paste0(data_dir, "/clustertracker/node_to_leaves.csv"))
    clustertracker_descendent_data_tmp <- node_to_leaves %>%
        mutate(method = "clustertracker") %>%
        rename("label" = leaf) %>%
        tidyr::separate(introduction_node, remove = F, into = c("introduction_location", "node_idx"), sep = "_node_") %>%
        filter(introduction_location == outbreak_annotation) %>%
        left_join(full_metadata_tmp, by = c("label" = "strain")) %>%
        mutate(tip = label) %>%
        dplyr::select(tip, introduction_node, method, label, location)

    # Take only introduction closest to tips (in case of re-introduction)
    clustertracker_descendent_data <- clustertracker_descendent_data_tmp %>%
        group_by(tip) %>%
        slice_max(order_by = introduction_node, n = 1)

    ### Get clusters from hiv-trace output ###

    hivtrace_descendent_data <- rbind(cluster_results, singleton_results) %>% 
        mutate(location = case_when(
            sample %in% outbreak_sample_names ~ outbreak_annotation, T ~ "Other")) %>%
        mutate(tip = sample, label = sample, method = "hivtrace") %>%
        dplyr::select(tip, introduction_node, method, label, location)


    ### Combine all results ###

    full_results_descendent_data <- rbind(
        mugration_descendent_data,
        clustertracker_descendent_data,
        hivtrace_descendent_data
    ) %>% rename("tip_location" = location)

    # inspect results
    tmp <- full_results_descendent_data %>%
        group_by(method, introduction_node) %>%
        summarize(n_descendents = n(), n_focal_region = sum(tip_location == outbreak_annotation))
    print(tmp %>% group_by(method) %>% summarize(n_descendents = sum(n_descendents)))

    # save results
    write.table(
        x = full_results_descendent_data, 
        file = paste0(data_dir, "/collated_results_per_introduction.txt"),
        row.names = F)

}
