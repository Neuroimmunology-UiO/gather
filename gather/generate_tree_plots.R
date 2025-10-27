# generate_tree_plots.R
options(warn = -1)

# Auto-install required packages if missing
required_packages <- c("dowser", "ggtree", "dplyr", "ggrepel", "RColorBrewer")

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
file_path      <- args[1]
file_path_heavy <- args[2]
germline_path  <- args[3]
output_dir     <- args[4]
igphyml_exec   <- args[5]

# Read in data
ExampleMixedDb <- read.delim(file_path, stringsAsFactors = FALSE)
ExampleMixedDb <- resolveLightChains(ExampleMixedDb)
ExampleMixedDb_heavy <- read.delim(file_path_heavy, stringsAsFactors = FALSE)

# Build custom palette
c_call_color_df <- ExampleMixedDb_heavy[, c("c_call", "color")]
c_call_color_df <- c_call_color_df[!is.na(c_call_color_df$c_call) & !is.na(c_call_color_df$color), ]
c_call_color_df_unique <- c_call_color_df[!duplicated(c_call_color_df$c_call), ]

custom_palette <- setNames(
  as.character(c_call_color_df_unique$color),
  as.character(c_call_color_df_unique$c_call)
)

# Make sure all c_call values present in the main DB are covered
all_calls <- unique(na.omit(ExampleMixedDb$c_call))
missing_calls <- setdiff(all_calls, names(custom_palette))
if (length(missing_calls) > 0) {
  extra_colors <- RColorBrewer::brewer.pal(
    max(3, length(missing_calls)), "Set3"
  )[seq_along(missing_calls)]
  custom_palette <- c(custom_palette, setNames(extra_colors, missing_calls))
}

# Load references
references <- readIMGT(germline_path)

ExampleMixedDb <- createGermlines(
  ExampleMixedDb,
  references = references,
  clone = "clone_subgroup_id",
  nproc = 1
)

# Format clones
clones <- formatClones(
  ExampleMixedDb,
  chain = "HL",
  nproc = 1,
  collapse = TRUE,
  add_count = TRUE,
  split_light = TRUE,
  minseq = 2,
  text_fields = "c_call",
  traits = "cell_id"
)

# Get trees
clones <- getTrees(
  clones,
  build = "igphyml",
  nproc = 1,
  partition = "hl",
  exec = igphyml_exec
)

# Scale branches
ExampleClones_m <- scaleBranches(clones, edge_type = "mutations")

# Plot trees (wrap to quiet upstream deprecation chatter from dependencies)
plots <- suppressWarnings(
  plotTrees(
    ExampleClones_m,
    scale = 10,
    tips = "c_call",
    tipsize = "collapse_count",
    palette = custom_palette,
    layout = "circular"
  )
)

# Add legends, theme, and cell_id labels
# NOTE: use tidy-eval with subset inside aes() to avoid aes_()/filter deprecation
treeplots <- lapply(plots, function(x) {
  x +
    ggtree::geom_text2(
      aes(label = cell_id, subset = isTip),
      hjust = -0.2, vjust = 0.5, size = 3
    ) +
    guides(
      color = guide_legend(title = "Cell types"),
      size  = guide_legend(title = "Collapsed cells")
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
        filename = file.path(output_dir, paste0("treeplot", i, ".pdf")),
        plot = treeplots[[i]],
        width = 7,
        height = 5
      )
    }
  }
}

save_plots(treeplots, save_all = TRUE)