# Load Libraries
library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)

setwd("ADD YOUR PATH")

out_dir <- "tree_plots_all"

# Create folder only if it does NOT exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# List Data
gene_of_interest <- c("EF1", "ITS2")

# =================== For Loop ===================
for (gene in gene_of_interest) {
  
  gene_treefile <- paste0("file_tree_all_sample/", gene, "_all_tree.treefile")
  treefile <- read.tree(gene_treefile)
  
  p <- ggtree(treefile)
  
  data <- p$data
  
  tip_labs <- data %>%
    filter(isTip == TRUE) %>%
    mutate(Clean_Labels = gsub("_", " ", label))
  
  # Drop Cannabis_sativa if present
  if ("Cannabis_sativa" %in% treefile$tip.label) {
    tree_sativa <- drop.tip(treefile, tip = "Cannabis_sativa")
  } else {
    tree_sativa <- treefile
  }
  
  tree_no_sativa <- ggtree(tree_sativa)
  
  # Find Fusarium tips
  fusarium_tips <- tree_sativa$tip.label[grepl("Fusarium", tree_sativa$tip.label)]
  
  # Define highlight node only if there are at least 2 Fusarium tips
  highlight_node <- NULL
  if (length(fusarium_tips) >= 2) {
    highlight_node <- MRCA(tree_sativa, fusarium_tips)
  }
  
  # Compute tip count and dynamic label size
  n_tips <- length(treefile$tip.label)
  label_size <- ifelse(n_tips > 150, 1.8,
                       ifelse(n_tips > 100, 2.2,
                              ifelse(n_tips > 60, 3, 4)))
  
  # Plot tree
  tree_microbes <- p %<+% tip_labs +
    geom_tiplab(
      aes(label = Clean_Labels),
      align = TRUE,
      linetype = 0,
      size = label_size,
      fontface = "italic"
    ) +
    ggtitle(paste0(gene, " Tree")) +
    theme(
      plot.title = element_text(size = 15, face = "bold")
    )
  
  # Add highlight only if highlight_node exists
  if (!is.null(highlight_node)) {
    tree_microbes <- tree_microbes +
      geom_hilight(node = highlight_node, fill = "lightblue", extend = 1.5)
  }
  
  save_name <- file.path(out_dir, paste0(gene, "_", sample_name, "_Tree.pdf"))
  ggsave(save_name, plot = tree_microbes, width = 40, height = 50, units = "cm", dpi = 300)
}
