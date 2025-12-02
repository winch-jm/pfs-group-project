# app.R
library(shiny)
library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(plotly)
library(tidyr)

#----------------- load data & build graph (once) -----------------#
color_palette <- function(n) Polychrome::palette36.colors(n)


edges <- read_csv("../data/graph_edges_mcf7.csv",
                  col_types = cols(
                    src    = col_character(),
                    dist   = col_character(),
                    weight = col_double()
                  ))

# igraph uses 'name' as vertex ID by default; we'll treat src/dst as labels
g <- graph_from_data_frame(edges, directed = FALSE)
E(g)$weight <- edges$weight

metadata   <- read.table("../data/siginfo_beta.txt",    sep = "\t", header = TRUE)
cp_info    <- read.table("../data/compoundinfo_beta.txt",
                         sep = "\t", header = TRUE,
                         fill = TRUE, quote = "", comment.char = "")
atc_info   <- read_csv("../data/drug2atc.csv")
umap_coords <- read_csv("../data/MCF7_umap.csv")

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

partition <- readr::read_csv("../data/partition_mcf7_fewer.csv",
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


# -------------------------------------------------------------------------
# UI
# -------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Drug Association Network Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput(
        "color_by",
        "Color nodes by:",
        choices = c(
          "community",
          "cell_iname",
          "tas",
          "atc_level1" = "level1",
          "atc_level2" = "level2",
          "atc_level3" = "level3",
          "atc_level4" = "level4",
          "moa",
          "pert_time",
          "pert_dose"
        ),
        selected = "community"
      ),
      # --- NEW FILTERS ---
      selectizeInput(
        "cell_filter",
        "Cell line:",
        choices  = sort(unique(drug_nodes$cell_iname)),
        selected = sort(unique(drug_nodes$cell_iname)),
        multiple = TRUE
      ),
      selectizeInput(
        "time_filter",
        "Perturbation time (h):",
        choices  = sort(unique(drug_nodes$pert_time)),
        selected = sort(unique(drug_nodes$pert_time)),
        multiple = TRUE
      ),
      selectizeInput(
        "dose_filter",
        "Perturbation dose:",
        choices  = sort(unique(drug_nodes$pert_dose)),
        selected = sort(unique(drug_nodes$pert_dose)),
        multiple = TRUE
      ),
      checkboxInput("drop_no_atc", "Drop nodes without ATC label", value = FALSE),
      sliderInput("edge_alpha", "Edge transparency:",
                  min = 0.0, max = 1.0, value = 0.7, step = 0.1),
      sliderInput("node_size", "Base node size:",
                  min = 1, max = 10, value = 5, step = 1),
      checkboxInput("show_edges", "Show drug-level edges", value = TRUE),
      selectInput(
        "atc_level",
        "ATC level:",
        choices = c("Level 1" = "level1",
                    "Level 2" = "level2",
                    "Level 3" = "level3",
                    "Level 4" = "level4"),
        selected = "level1"
      )
    ),
    
    mainPanel(
      fluidRow(
        column(
          width = 12,
          h4("Graph layout (drug nodes)"),
          plotlyOutput("graph_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("UMAP (drug nodes)"),
          plotlyOutput("umap_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 6,
          h4("Community sizes (drug-level)"),
          plotOutput("comm_bar", height = "250px")
        ),
        column(
          width = 6,
          h4("TAS distribution per community (drug-level)"),
          plotOutput("tas_distribution", height = "250px")
        )
      ),
      br(),
      # Row 2b: MoA vs TAS
      fluidRow(
        column(
          width = 12,
          h4("TAS by MoA (drug-level)"),
          plotOutput("moa_tas_plot", height = "500px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12,
          h4("ATC Enrichment"),
          plotOutput("atc_enrichment", height = "600px")
        )
      ),
      fluidRow(
        column(
          width = 12,
          h4("MoA Enrichment"),
          plotOutput("moa_enrichment", height = "600px")
        )
      )
    )
  )
)

# -------------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------------

server <- function(input, output, session) {
  # returns colors for categorical or numeric variables
  get_colors <- function(values) {
    if (is.numeric(values)) {
      # continuous variable (TAS)
      return(~values)
    } else {
      # categorical variable
      lvls <- unique(values)
      pal  <- color_palette(length(lvls))
      names(pal) <- lvls
      return(pal[values])
    }
  }
  
  # subset drug_nodes by ATC presence if requested
  filtered_nodes <- reactive({
    df <- drug_nodes
    
    # ATC filter
    if (input$drop_no_atc) {
      df <- df %>%
        filter(
          !is.na(level1) | !is.na(level2) |
            !is.na(level3) | !is.na(level4)
        )
    }
    
    # cell line filter
    if (!is.null(input$cell_filter) && length(input$cell_filter) > 0) {
      df <- df %>% filter(cell_iname %in% input$cell_filter)
    }
    
    # perturbation time filter
    if (!is.null(input$time_filter) && length(input$time_filter) > 0) {
      df <- df %>% filter(pert_time %in% input$time_filter)
    }
    
    # dose filter
    if (!is.null(input$dose_filter) && length(input$dose_filter) > 0) {
      df <- df %>% filter(pert_dose %in% input$dose_filter)
    }
    
    df
  })
  
  atc_data <- reactive({
    lvl <- input$atc_level   # "level1"..."level4"
    
    df <- filtered_nodes() %>%
      filter(!is.na(.data[[lvl]])) %>%
      select(community, atc = all_of(lvl))
    
    N <- nrow(df)
    
    enr <- df %>%
      count(community, atc, name = "k") %>%
      left_join(df %>% count(atc,       name = "K"), by = "atc") %>%
      left_join(df %>% count(community, name = "n"), by = "community") %>%
      mutate(
        p         = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
        score     = -log10(p + 1e-12),
        score_sig = ifelse(p < 0.05, score, NA_real_),
        score_clean = ifelse(is.na(score_sig), 0, score_sig)  # use 0 for “no enrichment”
      )
    
    # ---- order communities by *their* max enrichment (descending) ----
    row_ord <- enr %>%
      group_by(community) %>%
      summarise(row_max = max(score_clean, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(row_max)) %>%
      pull(community)
    
    # ---- order ATC classes by *their* max enrichment (descending) ----
    col_ord <- enr %>%
      group_by(atc) %>%
      summarise(col_max = max(score_clean, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(col_max)) %>%
      pull(atc)
    
    enr$community <- factor(enr$community, levels = row_ord)
    enr$atc       <- factor(enr$atc,       levels = col_ord)
    
    enr
  })
  # ---- UMAP plot: NO linked brushing, just interactive scatter ----
  output$umap_plot <- renderPlotly({
    df <- filtered_nodes()
    color_col <- input$color_by
    
    # --- size: make it clearly depend on n_sigs ---
    # e.g. log scaling so 1 vs 40 reps are visibly different
    df <- df %>%
      mutate(
        size_raw = log2(n_sigs + 1),             # 1 rep -> 1, 3 reps -> 2, etc.
        # normalize to a reasonable range, then multiply by base node_size
        point_size = input$node_size * (size_raw / max(size_raw, na.rm = TRUE)),
        point_size = ifelse(is.finite(point_size), point_size, input$node_size)
      )
    
    # hover text
    df <- df %>%
      mutate(
        hover = if (is.numeric(.data[[color_col]])) {
          sprintf(
            "%s<br>%s: %.3f<br>#sigs: %d",
            cmap_name, color_col, .data[[color_col]], n_sigs
          )
        } else {
          sprintf(
            "%s<br>%s: %s<br>#sigs: %d",
            cmap_name, color_col, as.character(.data[[color_col]]), n_sigs
          )
        }
      )
    
    # sentinel for unlabeled if you still want that behavior
    df <- df %>%
      mutate(
        color_mapped = ifelse(
          is.na(.data[[color_col]]) | .data[[color_col]] == "",
          "UNLABELED__GRAY__",
          as.character(.data[[color_col]])
        )
      )
    
    plot_ly(
      data   = df,
      x      = ~UMAP_1,
      y      = ~UMAP_2,
      type   = "scattergl",
      mode   = "markers",
      color  = ~color_mapped,
      text   = ~hover,
      hoverinfo = "text",
      # >>> THIS is the important part: size is mapped here, not inside marker$size
      size   = ~point_size,
      marker = list(opacity = 0.8),
      sizes  = c(4, 18)     # overall min/max in screen units
    ) %>%
      layout(
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2")
      )
  })
  
  
  # ---- Graph plot: drug nodes + drug-level edges ----
  output$graph_plot <- renderPlotly({
    df <- filtered_nodes()
    color_col <- input$color_by
    
    # size scaling identical to UMAP for consistency
    df <- df %>%
      mutate(
        size_raw = log2(n_sigs + 1),
        point_size = input$node_size * (size_raw / max(size_raw, na.rm = TRUE)),
        point_size = ifelse(is.finite(point_size), point_size, input$node_size),
        hover = if (is.numeric(.data[[color_col]])) {
          sprintf(
            "%s<br>%s: %.3f<br>#sigs: %d",
            cmap_name, color_col, .data[[color_col]], n_sigs
          )
        } else {
          sprintf(
            "%s<br>%s: %s<br>#sigs: %d",
            cmap_name, color_col, as.character(.data[[color_col]]), n_sigs
          )
        }
      )
    
    p <- plot_ly(type = "scattergl", mode = "markers")
    
    # edges (unchanged)
    if (isTRUE(input$show_edges)) {
      edges_df <- drug_edges %>%
        filter(src_agg %in% df$agg_id, dst_agg %in% df$agg_id) %>%
        left_join(df %>% select(agg_id, x0 = graph_x, y0 = graph_y),
                  by = c("src_agg" = "agg_id")) %>%
        left_join(df %>% select(agg_id, x1 = graph_x, y1 = graph_y),
                  by = c("dst_agg" = "agg_id"))
      
      if (nrow(edges_df) > 0) {
        p <- p %>%
          add_segments(
            data = edges_df,
            x    = ~x0,
            y    = ~y0,
            xend = ~x1,
            yend = ~y1,
            line = list(width = 0.3,
                        color = paste0("rgba(150,150,150,", input$edge_alpha, ")")),
            hoverinfo = "none",
            showlegend = FALSE
          )
      }
    }
    
    # drug nodes – note size argument!
    p <- p %>%
      add_markers(
        data  = df,
        x     = ~graph_x,
        y     = ~graph_y,
        mode  = "markers",
        color = ~.data[[color_col]],
        text  = ~hover,
        hoverinfo = "text",
        size  = ~point_size,
        marker = list(opacity = 0.9),
        sizes = c(4, 18),
        showlegend = TRUE
      ) %>%
      layout(
        xaxis = list(title = NULL, showticklabels = FALSE, zeroline = FALSE),
        yaxis = list(title = NULL, showticklabels = FALSE, zeroline = FALSE)
      )
    
    p
  })
  
  # ---- community size bar ----
  output$comm_bar <- renderPlot({
    filtered_nodes() %>%
      count(community) %>%
      ggplot(aes(x = community, y = n)) +
      geom_col() +
      labs(title = "Community Sizes (drug-level)",
           x = "Community", y = "# Drugs") +
      theme_minimal()
  })
  
  # ---- TAS distribution ----
  output$tas_distribution <- renderPlot({
    ggplot(filtered_nodes(), aes(x = community, y = tas, fill = community)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "TAS Distribution per Community (drug-level)",
        x = "Community",
        y = "Mean TAS per drug"
      ) +
      theme(legend.position = "none")
  })
  
  # ---- ATC enrichment heatmap ----
  output$atc_enrichment <- renderPlot({
    enr <- atc_data()
    
    level_label <- switch(
      input$atc_level,
      "level1" = "ATC Level 1",
      "level2" = "ATC Level 2",
      "level3" = "ATC Level 3",
      "level4" = "ATC Level 4"
    )
    
    ggplot(enr, aes(x = community, y = atc, fill = score_sig)) +
      geom_tile() +
      scale_fill_viridis_c(
        option   = "C",
        na.value = "grey95",
        name     = "-log10 p"
      ) +
      labs(
        title = paste("ATC Enrichment by Community:", level_label),
        x = "Community",
        y = level_label
      ) +
      theme_minimal()
  })
  
  # ---- MoA enrichment dotplot ----
  # ---- MoA enrichment with diagonalized ordering ----
  moa_data <- reactive({
    df <- filtered_nodes() %>%
      filter(!is.na(moa), moa != "") %>%
      select(community, moa)
    
    N <- nrow(df)
    
    enr <- df %>%
      count(community, moa, name = "k") %>%
      left_join(df %>% count(moa,       name = "K"), by = "moa") %>%
      left_join(df %>% count(community, name = "n"), by = "community") %>%
      mutate(
        p     = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
        score = -log10(p + 1e-12),
        fold  = (k / n) / (K / N)
      ) %>%
      filter(p < 0.05) %>%                     # keep significant MoAs
      group_by(community) %>%
      slice_max(order_by = score, n = 3, with_ties = FALSE) %>%  # top 3 per community
      ungroup()
    
    ## ---- 1. Order communities by their max MoA enrichment (increasing) ----
    comm_order_df <- enr %>%
      group_by(community) %>%
      summarize(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
      arrange(max_score) %>%                   # increasing; use desc() for decreasing
      mutate(comm_rank = row_number())
    
    comm_order <- comm_order_df$community
    
    ## ---- 2. Order MoAs by the community where they reach their max ----
    moa_order <- enr %>%
      group_by(moa) %>%
      summarize(
        max_comm      = community[which.max(score)],
        max_score_moa = max(score, na.rm = TRUE),
        .groups       = "drop"
      ) %>%
      mutate(
        comm_rank = match(max_comm, comm_order)  # position of that community
      ) %>%
      arrange(comm_rank, desc(max_score_moa)) %>% # so big blobs march along the diagonal
      pull(moa)
    
    enr$community <- factor(enr$community, levels = comm_order)
    enr$moa       <- factor(enr$moa,       levels = moa_order)
    
    enr
  })
  
  
  output$moa_enrichment <- renderPlot({
    enr <- moa_data()
    
    ggplot(enr, aes(x = community, y = moa)) +
      geom_point(
        aes(size = fold, color = score),
        alpha = 0.9
      ) +
      scale_color_viridis_c(
        option = "C",
        name   = "-log10 p"
      ) +
      scale_size(
        range = c(2, 8),
        name  = "Fold enrichment"
      ) +
      labs(
        title = "MoA Enrichment by Community (dotplot, p < 0.05)",
        x = "Community",
        y = "Mechanism of Action"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
  output$moa_tas_plot <- renderPlot({
    # use aggregated drug-level nodes, respecting ATC filter
    df <- filtered_nodes() %>%
      dplyr::filter(!is.na(moa), moa != "") 
    
    # keep MoAs with at least N drugs (tune this)
    min_drugs <- 5
    df <- df %>%
      add_count(moa, name = "n_drugs") %>%
      dplyr::filter(n_drugs >= min_drugs)
    
    # order MoAs by median TAS (high to low)
    df <- df %>%
      group_by(moa) %>%
      mutate(median_tas = median(tas, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(moa = reorder(moa, median_tas))
    
    ggplot(df, aes(x = tas, y = moa)) +
      geom_boxplot(
        aes(group = moa),
        fill = "grey90",
        color = "grey40",
        outlier.size = 0.6
      ) +
      geom_jitter(
        width = 0,
        height = 0.1,
        alpha = 0.4,
        size  = 0.8,
        aes(color = community)
      ) +
      scale_color_discrete(name = "Community") +
      labs(
        title = "TAS by Mechanism of Action (drug-level)",
        x = "Mean TAS per drug",
        y = "Mechanism of Action"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 7),
        legend.position = "right"
      )
  })
  
  
}

shinyApp(ui, server)

