# app.R
library(shiny)
library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(plotly)
library(tidyr)
library(viridisLite)
library(ggridges)
library(ggrepel)
library(ggalluvial)

#------------------------- data paths -----------------------------#

graphEdgeFile   <- "../test/gridsearch/pfsDemo_graph_edges_k15.csv"
partitionFile   <- "../test/gridsearch/pfsDemo_partition_k15_g50.csv"
sigInfoFile     <- "../data/siginfo_beta.txt"
cpInfoFile      <- "../data/compoundinfo_beta.txt"
atcInfoFile     <- "../data/drug2atc.csv"
umapCoordsFile  <- "../data/MCF7_umap.csv"



#----------------- load data & build graph (once) -----------------#
color_palette <- function(n) Polychrome::palette36.colors(n)


edges <- read_csv(graphEdgeFile,
                  col_types = cols(
                    src    = col_character(),
                    dist   = col_character(),
                    weight = col_double()
                  ))

# igraph uses 'name' as vertex ID by default; we'll treat src/dst as labels
g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$weight

metadata   <- read.table(sigInfoFile,    sep = "\t", header = TRUE)
cp_info    <- read.table(cpInfoFile,
                         sep = "\t", header = TRUE,
                         fill = TRUE, quote = "", comment.char = "")
atc_info   <- read_csv(atcInfoFile)
umap_coords <- read_csv(umapCoordsFile)

# --- de-duplicate cp_info & atc_info to avoid many-to-many joins --- #
cp_info_unique <- cp_info %>%
  group_by(cmap_name, pert_id) %>%
  slice(1) %>%    # keep first row per (cmap_name, pert_id)
  ungroup()

atc_info_unique <- atc_info %>%
  group_by(inchi_key) %>%
  slice(1) %>%    # keep one ATC row per compound
  ungroup()

# --- join metadata with compound & umap & ATC --- #
metadata <- metadata %>%
  left_join(cp_info_unique,  by = c("cmap_name", "pert_id")) %>%
  left_join(umap_coords,     by = "sig_id") %>%
  left_join(atc_info_unique, by = "inchi_key")

partition <- readr::read_csv(partitionFile,
                             col_types = cols(
                               node      = col_character(),
                               node_id   = col_character(),
                               community = col_integer()
                             ))

# make community a factor
partition <- partition %>%
  mutate(
    node      = as.character(node),
    node_id   = as.character(node_id),
    community = factor(community)
  )

# Merge partition + metadata
node_info <- partition %>%
  left_join(metadata, by = c("node_id" = "sig_id")) %>%
  mutate(
    atc_class = level1,
    atc_desc  = level1_description
  )

print(unique(node_info$moa))

# reorder communities by size (signature-level, mostly just for consistency)
comm_order <- node_info %>%
  count(community, name = "n") %>%
  arrange(desc(n)) %>%
  pull(community)

node_info$community <- factor(node_info$community, levels = comm_order)

# attach attributes to original igraph vertices (optional now)
V(g)$community  <- node_info$community[match(V(g)$name, node_info$node)]
V(g)$cell_iname <- node_info$cell_iname[match(V(g)$name, node_info$node)]
V(g)$pert_id    <- node_info$pert_id[match(V(g)$name, node_info$node)]
V(g)$tas        <- node_info$tas[match(V(g)$name, node_info$node)]
V(g)$atc_class  <- node_info$level1[match(V(g)$name, node_info$node)]
V(g)$atc_desc   <- node_info$level1_description[match(V(g)$name, node_info$node)]
V(g)$moa        <- node_info$moa[match(V(g)$name, node_info$node)]
V(g)$drug       <- node_info$cmap_name[match(V(g)$name, node_info$node)]

cat("ATC class distribution (vertex-level):\n")
print(table(V(g)$atc_class))
cat("vcount(g):", vcount(g), "\n")

# layout on original graph (signature-level) – used to derive drug-level positions
graph_layout <- create_layout(g, layout = "fr")

layout_df <- graph_layout %>%
  as.data.frame() %>%
  rename(graph_x = x, graph_y = y) %>%
  select(name, graph_x, graph_y)

# attach layout coords to node_info
node_info <- node_info %>%
  left_join(layout_df, by = c("node" = "name"))

# -------- aggregate replicates to drug × community level -------- #
drug_nodes <- node_info %>%
  group_by(cmap_name, community) %>%
  summarise(
    node_ids   = paste(node, collapse = ";"),
    node_id    = first(node_id),
    # average coords
    UMAP_1     = mean(UMAP_1, na.rm = TRUE),
    UMAP_2     = mean(UMAP_2, na.rm = TRUE),
    graph_x    = mean(graph_x, na.rm = TRUE),
    graph_y    = mean(graph_y, na.rm = TRUE),
    # TAS and replicate count
    tas        = mean(tas, na.rm = TRUE),
    n_sigs     = n(),
    # annotations
    cell_iname = first(cell_iname),
    pert_time  = first(pert_time),     # <-- NEW
    pert_dose  = first(pert_dose),     # <-- NEW
    level1     = first(level1),
    level2     = first(level2),
    level3     = first(level3),
    level4     = first(level4),
    atc_desc   = first(level1_description),
    moa        = first(moa),
    .groups    = "drop"
  ) %>%
  mutate(
    atc_class = level1,
    agg_id    = row_number()  # unique ID for aggregated nodes
  )

# reorder communities by drug-level size
comm_order_drug <- drug_nodes %>%
  count(community, name = "n") %>%
  arrange(desc(n)) %>%
  pull(community)

drug_nodes$community <- factor(drug_nodes$community, levels = comm_order_drug)

cat("signature-level layout rows:", nrow(graph_layout), "\n")
cat("NA coords:", sum(is.na(graph_layout$x) | is.na(graph_layout$y)), "\n")

# ---- build drug-level edge list from signature-level edges ---- #

# map each original node to an agg_id (drug × community)
node_to_agg <- node_info %>%
  select(node, cmap_name, community) %>%
  inner_join(drug_nodes %>% select(cmap_name, community, agg_id),
             by = c("cmap_name", "community"))

drug_edges <- edges %>%
  rename(src_node = src, dst_node = dist) %>%
  left_join(node_to_agg, by = c("src_node" = "node")) %>%
  rename(src_agg = agg_id) %>%
  left_join(node_to_agg, by = c("dst_node" = "node")) %>%
  rename(dst_agg = agg_id) %>%
  filter(!is.na(src_agg), !is.na(dst_agg), src_agg != dst_agg) %>%
  mutate(
    a = pmin(src_agg, dst_agg),
    b = pmax(src_agg, dst_agg)
  ) %>%
  group_by(a, b) %>%
  summarise(weight = mean(weight), .groups = "drop") %>%
  rename(src_agg = a, dst_agg = b)

cat("drug-level edges:", nrow(drug_edges), "\n")


# ---- signature-level edges annotated with community on each endpoint ----
edge_comm <- edges %>%
  # join community for source node
  left_join(
    partition %>% select(node, community),
    by = c("src" = "node")
  ) %>%
  rename(comm_src = community) %>%
  # join community for destination node
  left_join(
    partition %>% select(node, community),
    by = c("dist" = "node")
  ) %>%
  rename(comm_dst = community) %>%
  mutate(
    comm_src = as.character(comm_src),
    comm_dst = as.character(comm_dst)
  )