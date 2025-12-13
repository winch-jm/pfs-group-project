library(tidyverse)
library(igraph)

#-----------------------------
# Paths & metadata
#-----------------------------
grid_dir   <- "../test/gridsearch/"
siginfo    <- "../data/siginfo_beta.txt"

metadata   <- read.table("../data/siginfo_beta.txt",    sep = "\t", header = TRUE)
cp_info    <- read.table("../data/compoundinfo_beta.txt",
                         sep = "\t", header = TRUE,
                         fill = TRUE, quote = "", comment.char = "")
atc_info   <- read_csv("../data/drug2atc.csv")
umap_coords <- read_csv("../data/MCF7_umap.csv")

# --- de-duplicate cp_info & atc_info to avoid many-to-many joins --- #
cp_info_unique <- cp_info %>%
  group_by(cmap_name, pert_id) %>%
  slice(1) %>%
  ungroup()

atc_info_unique <- atc_info %>%
  group_by(inchi_key) %>%
  slice(1) %>%
  ungroup()

# --- join metadata with compound & umap & ATC --- #
metadata <- metadata %>%
  left_join(cp_info_unique,  by = c("cmap_name", "pert_id")) %>%
  left_join(umap_coords,     by = "sig_id") %>%
  left_join(atc_info_unique, by = "inchi_key")

#-----------------------------
# helper: compute MoA purity (entropy-based)
#-----------------------------
compute_moa_purity <- function(partition_df, metadata_df) {
  # partition_df: columns node, node_id, community
  # metadata_df:  columns sig_id, moa
  
  comm_moa <- partition_df %>%
    left_join(metadata_df, by = c("node_id" = "sig_id")) %>%
    filter(!is.na(moa), moa != "") %>%        # ignore unlabeled
    group_by(community, moa) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(community) %>%
    mutate(
      total = sum(n),
      p     = n / total
    ) %>%
    summarise(
      H     = -sum(p * log(p)),   # Shannon entropy (natural log)
      K     = n(),                # # distinct MoAs in community
      n_lab = first(total),       # # labeled members in this community
      .groups = "drop"
    ) %>%
    mutate(
      H_max  = ifelse(K > 1, log(K), 0),
      purity = case_when(
        H_max > 0 ~ 1 - H / H_max,   # normalized entropy
        K == 1    ~ 1,               # pure single-MoA community
        TRUE      ~ NA_real_
      )
    )
  
  # Weighted average purity across communities, weighted by # labeled members
  comm_moa %>%
    filter(!is.na(purity), n_lab > 0) %>%
    summarise(
      moa_purity = sum(purity * n_lab) / sum(n_lab)
    ) %>%
    pull(moa_purity)
}

#-----------------------------
# Get all partition files
#-----------------------------
part_files <- list.files(
  grid_dir,
  pattern = "^partition_mcf7_k\\d+_g\\d+\\.csv$",
  full.names = TRUE
)

# Pre-load all graphs keyed by k (10,15,20,25)
k_values <- c(10, 15, 20, 25)

graphs_by_k <- map(
  k_values,
  ~ {
    edge_path <- file.path(grid_dir, sprintf("pfsDemo_graph_edges_k%d.csv", .x))
    edges <- readr::read_csv(
      edge_path,
      col_types = cols(
        src    = col_character(),
        dist   = col_character(),
        weight = col_double()
      )
    )
    g <- graph_from_data_frame(edges, directed = FALSE)
    E(g)$weight <- edges$weight
    g
  }
) %>%
  set_names(as.character(k_values))

#-----------------------------
# Loop over all (k, gamma) partitions
#-----------------------------
run_stats <- map_dfr(part_files, function(path) {
  fname <- basename(path)
  
  # parse k and gamma*100 from filename
  m <- stringr::str_match(fname, "partition_mcf7_k(\\d+)_g(\\d+)\\.csv")
  k      <- as.integer(m[2])
  g100   <- as.integer(m[3])
  gamma  <- g100 / 100   # original γ values: 0.25, 0.50, ...
  
  # read partition
  part <- readr::read_csv(
    path,
    col_types = cols(
      node      = col_character(),
      node_id   = col_character(),
      community = col_integer()
    )
  )
  
  # get graph for this k
  g <- graphs_by_k[[as.character(k)]]
  
  # map nodes to membership vector in the order of V(g)$name
  memb <- part$community[match(V(g)$name, part$node)]
  valid <- !is.na(memb)
  g_sub <- induced_subgraph(g, vids = V(g)[valid])
  memb_sub <- memb[valid]
  
  # modularity
  Q <- modularity(g_sub, membership = memb_sub, weights = E(g_sub)$weight)
  
  # MoA purity
  moa_purity <- compute_moa_purity(part, metadata)
  
  tibble(
    k          = k,
    gamma      = gamma,
    modularity = Q,
    moa_purity = moa_purity
  )
})
print(max(run_stats$modularity))
#-----------------------------
# Normalize modularity & build composite
#-----------------------------
mod_range <- range(run_stats$modularity, na.rm = TRUE)
grid_stats <- run_stats %>%
  mutate(
    norm_modularity = (modularity - mod_range[1]) /
      (mod_range[2] - mod_range[1]),
    composite_score = 0.5 * (moa_purity + norm_modularity),
    
    # ordered factor for gamma so y-ticks show original values
    gamma_f = factor(
      gamma,
      levels = sort(unique(gamma)),   # numeric order
      labels = sort(unique(gamma))    # printed labels (0.25, 0.50, ...)
    )
  )

print(grid_stats)

#-----------------------------
# Plot: heatmap over (k, gamma)
#-----------------------------
ggplot(grid_stats, aes(x = factor(k), y = gamma_f, fill = composite_score)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Composite\n(MoA purity + norm Q)/2") +
  labs(
    title = "Grid search over k and resolution",
    subtitle = "Composite of MoA entropy-purity and normalized modularity",
    x = "k (KNN neighbors)",
    y = "Resolution γ"
  ) +
  theme_minimal(base_size = 13)
