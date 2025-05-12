# generate_tree_plots.R

# Auto-install required packages if missing
required_packages <- c("dowser", "ggtree", "dplyr", "ggrepel")

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

install_if_missing(required_packages)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]
file_path_heavy <- args[2]
germline_path <- args[3]
output_dir <- args[4]
igphyml_exec <- args[5]

# Read in data
ExampleMixedDb <- read.delim(file_path, stringsAsFactors = FALSE)
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)

ExampleMixedDb_heavy <- read.delim(file_path_heavy, stringsAsFactors = FALSE)
c_call_color_df <- ExampleMixedDb_heavy[, c("c_call", "color")]
c_call_color_df_unique <- c_call_color_df[!duplicated(c_call_color_df$c_call), ]
custom_palette <- setNames(c_call_color_df_unique$color, c_call_color_df_unique$c_call)

# Load references
references <- readIMGT(germline_path)
ExampleMixedDb <- createGermlines(ExampleMixedDb, references = references, clone = "clone_subgroup_id", nproc = 1)

# Format clones
clones <- formatClones(ExampleMixedDb, chain = "HL", nproc = 1, collapse = TRUE, add_count = TRUE,
                       split_light = TRUE, minseq = 2, text_fields = "c_call", traits = "cell_id")

# Get trees
clones <- getTrees(clones, build = "igphyml", nproc = 1, partition = "hl", exec = igphyml_exec)

# Scale branches
ExampleClones_m <- scaleBranches(clones, edge_type = "mutations")

# Plot trees
plots <- plotTrees(ExampleClones_m, scale = 10, tips = "c_call",
                   tipsize = "collapse_count", palette = custom_palette, layout = "circular")

treeplots <- lapply(plots, function(x) {
  x + guides(
    color = guide_legend(title = "Cell types"),
    size = guide_legend(title = "Collapsed cells")
  ) +
    theme(
      text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
})

# Save plots
save_plots <- function(treeplots, save_all = TRUE) {
  if (save_all) {
    for (i in seq_along(treeplots)) {
      ggsave(
        filename = file.path(output_dir, paste0("treeplot", i, ".png")),
        plot = treeplots[[i]],
        dpi = 200
      )
    }
  }
}

save_plots(treeplots, save_all = TRUE)