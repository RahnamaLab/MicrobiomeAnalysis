# ============================================================
# Load Libraries
# ============================================================
library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)

setwd("ADD YOUR PATH")
out_dir <- "tree_plots_per_sample"

# Create folder only if it does NOT exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Genes to analyze
gene_of_interest <- c("EF1", "ITS2")

# ============================================================
# FUNCTION: Choose layout depending on gene
# ============================================================
choose_layout <- function(gene, tree) {
  if (gene == "ITS2") {
    return(ggtree(tree, layout = "circular"))
  } else {
    return(ggtree(tree))
  }
}

# ============================================================
# FUNCTION: dynamically scale tip-label size
# ============================================================
dynamic_label_size <- function(n_tips){
  if (n_tips > 300) return(0.75)
  if (n_tips > 150) return(1.8)
  if (n_tips > 100) return(2.2)
  if (n_tips > 60)  return(3)
  if (n_tips > 30)  return(4)
  return(5)
}

# ============================================================
# SAMPLE-LEVEL TREES
# ============================================================
tree_dir <- "file_tree_per_sample"

for (gene in gene_of_interest) {
  
  file_pattern <- paste0("^", gene, ".*\\.treefile$")
  
  tree_files <- list.files(
    path = tree_dir,
    pattern = file_pattern,
    full.names = TRUE
  )
  
  print(tree_files)
  
  for (f in tree_files) {
    
    cat("Processing:", f, "\n")
    
    sample_name <- basename(f)
    sample_name <- sub(paste0("^", gene, "_"), "", sample_name)
    sample_name <- sub("\\.treefile$", "", sample_name)
    
    tr <- read.tree(f)
    
    # Remove Cannabis_sativa tips if present
    cannabis_tips <- tr$tip.label[grepl("Cannabis_sativa", tr$tip.label)]
    if (length(cannabis_tips) > 0) {
      tr <- drop.tip(tr, cannabis_tips)
    }
    
    # Choose layout
    p <- choose_layout(gene, tr)
    
    data <- p$data
    tip_labs <- data %>%
      filter(isTip) %>%
      mutate(Clean_Labels = gsub("_", " ", label))
    
    n_tips <- length(tr$tip.label)
    label_size <- dynamic_label_size(n_tips)
    
    # Find Fusarium tips by tip labels, not node ids
    fusarium_tips <- tr$tip.label[grepl("Fusarium", tr$tip.label)]
    
    highlight_node <- NULL
    if (length(fusarium_tips) >= 2) {
      highlight_node <- MRCA(tr, fusarium_tips)
    }
    
    tree_microbes <- p %<+% tip_labs +
      geom_tiplab(
        aes(label = Clean_Labels),
        align = TRUE,
        linetype = 1,
        size = label_size,
        fontface = "italic"
      ) +
      ggtitle(paste0(gene, " ", sample_name, " Tree")) +
      theme(plot.title = element_text(size = 25, face = "bold"))
    
    if (!is.null(highlight_node)) {
      tree_microbes <- tree_microbes +
        geom_hilight(node = highlight_node, fill = "lightblue", extend = 1.5)
    }
    
    if (gene == "EF1") {
      tree_microbes <- tree_microbes + coord_cartesian(xlim = c(0, 2))
    }
    
    save_name <- file.path(out_dir, paste0(gene, "_", sample_name, "_Tree.pdf"))
    ggsave(
      save_name,
      plot = tree_microbes,
      width = 50,
      height = 50,
      units = "cm",
      dpi = 300
    )
  }
}
