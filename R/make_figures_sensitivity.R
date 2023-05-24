# This script is to make figures showing the results including sensitivity analyses

library(tidyverse)
library(ggplot2)
library(ggtree)
library(gridExtra)
library(ggpubr)

# Analysis-specific parameters
data_dirs <- list(
    "clean_results/confa_sensitivity/results_replicate_1", 
    "clean_results/confa_sensitivity/results_replicate_2", 
    "clean_results/confa_sensitivity/results_replicate_3", 
    "clean_results/urmc/results_2023-05-15T13:15:38.785201-07:00")
outbreak_annotations <- list(
    "CONF_A", 
    "CONF_A", 
    "CONF_A", 
    "URMC")
focal_region_annotation <- list("Massachusetts", "URMC")
focal_outbreak_samples <- list(
    "boston" = c(
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
        "MT520310.1"),
    urmc = c(
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
)

is_first <- T
for (i in 1:length(data_dirs)) {
    # Load collated results
    results_per_sample_tmp <- tryCatch(
        { read.delim(paste0(data_dirs[[i]], "/collated_results_per_sample_with_beast.txt"), sep = " ") }, 
        error = function(e) {
            print("Loading results with BEAST failed, looking for results without BEAST")
            return(read.delim(paste0(data_dirs[[i]], "/collated_results_per_sample.txt"), sep = " "))
        }
    )
    results_per_introduction_tmp <- tryCatch(
        { read.delim(paste0(data_dirs[[i]], "/collated_results_per_introduction_with_beast.txt"), sep = " ") }, 
        error = function(e) {
            print("Loading results with BEAST failed, looking for results without BEAST")
            return(read.delim(paste0(data_dirs[[i]], "/collated_results_per_introduction.txt"), sep = " "))
        }
    )
    
    results_per_sample_tmp <- results_per_sample_tmp %>% 
            mutate(
                analysis = i,
                outbreak = case_when(
                    outbreak_annotations[[i]] == "CONF_A" ~ "SARS-CoV-2 superspreading",
                    outbreak_annotations[[i]] == "URMC" ~ "K. aerogenes in CICU"),
                is_focal = sample %in% c(focal_outbreak_samples[[1]], focal_outbreak_samples[[2]]),  # this assumes annotations don't overlap
                method = factor(
                    method, 
                    levels = c("hivtrace", "clustertracker", "mugration", "beast"),
                    labels = c("HIV-\nTRACE", "ClusterTracker", "augur", "BEAST")))

    results_per_introduction_tmp <- results_per_introduction_tmp %>%
            mutate(
                analysis = i,
                outbreak = case_when(
                    outbreak_annotations[[i]] == "CONF_A" ~ "SARS-CoV-2 superspreading",
                    outbreak_annotations[[i]] == "URMC" ~ "K. aerogenes in CICU"),
                is_focal = label %in% c(focal_outbreak_samples[[1]], focal_outbreak_samples[[2]]),  # this assumes annotations don't overlap
                is_focal_region = tip_location %in% c(focal_region_annotation[[1]], focal_region_annotation[[2]]),  # this assumes annotations don't overlap
                method = factor(
                    method, 
                    levels = c("hivtrace", "clustertracker_no_asr", "clustertracker", "mugration", "beast"),
                    labels = c("HIV-\nTRACE", "ClusterTracker", "ClusterTracker\nASR", "augur", "BEAST")))

    if (is_first) {
        results_per_sample <- results_per_sample_tmp
        results_per_introduction <- results_per_introduction_tmp
        is_first <- F
    } else {
        results_per_sample <- rbind(results_per_sample, results_per_sample_tmp)
        results_per_introduction <- rbind(results_per_introduction, results_per_introduction_tmp)
    }
}


### Plot number of clusters by method ###

outbreak_clusters_by_method <- results_per_sample %>%
    filter(is_focal) %>%
    group_by(method, outbreak, analysis) %>%
    summarise(n_samples = n(), n_clusters = length(unique(introduction_node)))
outbreak_clusters_by_method_summary <- outbreak_clusters_by_method %>%
    group_by(method, outbreak) %>%
    summarize(
        mean_n_clusters = mean(n_clusters), 
        min_n_clusters = min(n_clusters),
        max_n_clusters = max(n_clusters)) %>%
    mutate(
        min_n_clusters_to_plot = case_when(
            outbreak == "K. aerogenes in CICU" ~ NA,
            method == "BEAST" ~ NA,
            T ~ min_n_clusters),
        max_n_clusters_to_plot = case_when(
            outbreak == "K. aerogenes in CICU" ~ NA,
            method == "BEAST" ~ NA,
            T ~ max_n_clusters))  # no replicates run, don't show errbars

outbreak_clusters_by_method_plot <- ggplot(data = outbreak_clusters_by_method_summary %>% filter(outbreak == "SARS-CoV-2 superspreading")) +
    geom_col(aes(x = method, y = mean_n_clusters)) +
    geom_errorbar(aes(x = method, ymin = min_n_clusters_to_plot, ymax = max_n_clusters_to_plot)) +
    # facet_grid(. ~ outbreak) +
    theme_bw() + 
    labs(x = element_blank(), y = "# clusters with outbreak samples") +
    geom_hline(aes(yintercept = 1), linetype = "dashed")

outbreak_clusters_by_method_plot


### Plot size of sampled outbreak ###

focal_outbreak_introductions <- results_per_introduction %>%
    filter(is_focal) %>%
    dplyr::select(introduction_node, method, outbreak, analysis) %>%
    unique()

outbreak_spread_by_method <- results_per_introduction %>%
    right_join(focal_outbreak_introductions) %>%  # only consider samples in focal outbreak cluster
    group_by(method, outbreak, analysis) %>%
    summarize(
        n_introduction_nodes = length(unique(introduction_node)),
        `Other` = sum(!is_focal_region & !is_focal),
        `Outbreak` = sum(is_focal),
        `Focal region` = sum(is_focal_region & !is_focal)) %>%
    pivot_longer(cols = c(`Other`, `Outbreak`, `Focal region`), names_to = "Sample type") %>%
    mutate(`Sample type` = factor(
        `Sample type`, 
        levels = c("Outbreak", "Focal region", "Other"),
        labels = c("Outbreak", "Focal region", "Other (divisional) region")))

outbreak_spread_by_method_summary <- outbreak_spread_by_method %>%
    group_by(method, outbreak, `Sample type`) %>%
    summarize(
        mean_value = mean(value),
        min_value = min(value),
        max_value = max(value)) %>%
    mutate(
        min_value_to_plot = case_when(
            outbreak == "K. aerogenes in CICU" ~ NA,
            method == "BEAST" ~ NA,
            T ~ min_value),
        max_value_to_plot = case_when(
            outbreak == "K. aerogenes in CICU" ~ NA,
            method == "BEAST" ~ NA,
            T ~ max_value)) # no replicates run, don't show errbars

shared_fill_scale <- scale_fill_manual(
    values = list(
        "Outbreak" = "#00BA38", 
        "Focal region" = "#619CFF", 
        "Other (divisional) region" = "#F8766D")) 


outbreak_spread_by_method_plot <- ggplot(
    data = outbreak_spread_by_method_summary %>% 
        filter(method != "ClusterTracker\nASR", outbreak == "SARS-CoV-2 superspreading")) +
    geom_bar(aes(x = method, y = mean_value, fill = `Sample type`), position = "dodge", stat = "identity") +
    geom_errorbar(aes(x = method, ymin = min_value_to_plot, ymax = max_value_to_plot), position = position_dodge2()) +
    theme_bw() +
    shared_fill_scale +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 6)) +
    # facet_grid(. ~ outbreak) +
    # theme(legend.position = c(0.335, 0.65), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25)) +
    labs(x = element_blank(), y = "Samples in the outbreak cluster(s)")

outbreak_spread_by_method_plot


### Plot all in one ###
legend_plot <- ggpubr::as_ggplot(ggpubr::get_legend(outbreak_spread_by_method_plot))
m <- matrix(NA, 2, 2)
m[1, 1] <- 1
m[1, 2] <- 2
m[2, 2] <- 3
png(
    filename = "clean_results/outbreak_comparison_sensitivity_sars_cov_2_only.png", 
    width = 6.5, height = 3.5, units = "in", res = 200)
gridExtra::grid.arrange(
    outbreak_clusters_by_method_plot,
    outbreak_spread_by_method_plot + theme(legend.position = "none"),
    legend_plot,
    layout_matrix = m,
    heights = c(1, 0.1))
dev.off()
