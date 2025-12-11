library(shiny)
library(ggalluvial)
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
  
  observeEvent(filtered_nodes(), {
    df <- filtered_nodes()
    comm_choices <- sort(unique(as.character(df$community)))
    
    # avoid empty choices if everything is filtered out
    if (length(comm_choices) == 0) {
      updateSelectInput(
        session, "selected_comm",
        choices = "",
        selected = ""
      )
    } else {
      updateSelectInput(
        session, "selected_comm",
        choices = comm_choices,
        selected = if (!is.null(input$selected_comm) &&
                       input$selected_comm %in% comm_choices) {
          input$selected_comm
        } else {
          comm_choices[1]
        }
      )
    }
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
    
    # size scaling (keep your existing logic)
    df <- df %>%
      mutate(
        size_raw = log2(n_sigs + 1)
      )
    
    # common hover text
    df <- df %>%
      mutate(
        hover = if (is.numeric(.data[[color_col]])) {
          sprintf("%s<br>%s: %.3f<br>#sigs: %d",
                  cmap_name, color_col, .data[[color_col]], n_sigs)
        } else {
          sprintf("%s<br>%s: %s<br>#sigs: %d",
                  cmap_name, color_col, as.character(.data[[color_col]]), n_sigs)
        }
      )
    
    ## ---------- TAS as continuous Viridis colorbar ----------
    if (color_col == "tas") {
      pal <- viridisLite::viridis(256)  # vector of hex colors
      
      return(
        plot_ly(
          data  = df,
          x     = ~UMAP_1,
          y     = ~UMAP_2,
          type  = "scattergl",
          mode  = "markers",
          text  = ~hover,
          hoverinfo = "text",
          color = ~tas,
          colors = pal,
          size   = ~size_raw,
          sizes  = c(input$node_size, input$node_size * 4),
          marker = list(
            opacity  = 0.8,
            colorbar = list(title = "TAS")
          )
        )
      )
    }
    
    ## ---------- everything else (categorical / other numeric) ----------
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
      size   = ~size_raw,
      sizes  = c(input$node_size, input$node_size * 4),
      marker = list(opacity = 0.8)
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
            line = list(
              width = 0.3,
              color = paste0("rgba(150,150,150,", input$edge_alpha, ")")
            ),
            hoverinfo = "none",
            showlegend = FALSE
          )
      }
    }
    
    # ---------- TAS: continuous colorbar ----------
    if (color_col == "tas") {
      p <- p %>%
        add_markers(
          data  = df,
          x     = ~graph_x,
          y     = ~graph_y,
          text  = ~hover,
          hoverinfo = "text",
          color = ~tas,
          colors = "Viridis",
          size   = ~size_raw,
          sizes  = c(input$node_size, input$node_size * 4),
          marker = list(
            opacity  = 0.9,
            colorbar = list(title = "TAS")
          ),
          showlegend = FALSE
        )
    } else {
      # ---------- default: categorical / discrete colors ----------
      p <- p %>%
        add_markers(
          data  = df,
          x     = ~graph_x,
          y     = ~graph_y,
          mode  = "markers",
          color = ~.data[[color_col]],
          text  = ~hover,
          hoverinfo = "text",
          size   = ~size_raw,
          sizes  = c(input$node_size, input$node_size * 4),
          marker = list(opacity = 0.9),
          showlegend = TRUE
        )
    }
    
    p %>%
      layout(
        xaxis = list(title = NULL, showticklabels = FALSE, zeroline = FALSE),
        yaxis = list(title = NULL, showticklabels = FALSE, zeroline = FALSE)
      )
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
    df <- filtered_nodes() %>%
      group_by(community) %>%
      mutate(mean_tas = mean(tas, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(community = reorder(community, -mean_tas))   # descending
    
    ggplot(df, aes(x = community, y = tas, fill = community)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "TAS Distribution per Community (drug-level)",
        x = "Community (sorted by mean TAS)",
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
  
  moa_entropy_data <- reactive({
    df <- filtered_nodes() %>%
      # keep only labeled MoAs
      dplyr::filter(!is.na(moa)) %>%
      dplyr::mutate(
        moa_clean = trimws(moa)
      ) %>%
      dplyr::filter(moa_clean != "\"\"")
    
    if (nrow(df) == 0) return(NULL)
    
    # counts of labeled MoAs per community
    comm_moa_counts <- df %>%
      dplyr::count(community, moa_clean, name = "n")
    
    # per-community entropy + dominant MoA
    comm_moa_counts %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        n_drugs = sum(n),
        k_moa   = dplyr::n(),
        # normalized Shannon entropy: 0 = pure, 1 = maximally mixed
        entropy = {
          p <- n / sum(n)
          H <- -sum(p * log(p))
          if (dplyr::n() > 1) H / log(dplyr::n()) else 0
        },
        dominant_moa = moa_clean[which.max(n)],
        dom_frac     = max(n / sum(n)),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_drugs > 0)
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
  
  output$tas_vs_size <- renderPlot({
    df <- filtered_nodes() %>%
      group_by(community) %>%
      summarise(
        n_drugs  = n(),
        mean_tas = mean(tas, na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      filter(!is.na(mean_tas))
    
    # mark only extremes
    df <- df %>%
      mutate(
        label = case_when(
          # mean_tas >= quantile(mean_tas, 0.90) ~ as.character(community),
          # mean_tas <= quantile(mean_tas, 0.10) ~ as.character(community),
          # n_drugs  >= quantile(n_drugs, 0.90)  ~ as.character(community),
          TRUE ~ ""
        )
      )
    
    ggplot(df, aes(x = n_drugs, y = mean_tas)) +
      geom_point(aes(color = mean_tas), size = 2.8, alpha = 0.9) +
      geom_smooth(
        method = "gam",
        formula = y ~ s(x, k = 5),
        se = FALSE,
        color = "black",
        linetype = "dashed"
      ) +
      geom_text(aes(label = label), vjust = -0.6, size = 3) +
      scale_x_log10(
        breaks = c(1, 5, 10, 25, 50, 100, 250, 500),
        labels = scales::comma
      ) +
      scale_color_viridis_c(name = "Mean TAS") +
      labs(
        title = "Mean TAS vs Community Size (drug-level)",
        x = "# Drugs in Community (log10 scale)",
        y = "Mean TAS"
      ) +
      theme_minimal(base_size = 13)
  })
  
  # which community was last clicked in UMAP or graph?
  selected_community <- reactive({
    req(input$selected_comm)
    input$selected_comm
  })
  
  
  
  community_data <- reactive({
    comm <- selected_community()
    req(comm)        # wait until something is clicked
    df <- filtered_nodes()
    df %>% dplyr::filter(community == comm)
  })
  
  output$community_summary <- renderTable({
    df <- community_data()
    
    tibble::tibble(
      Community          = as.character(df$community[1]),
      `# drugs`          = dplyr::n_distinct(df$cmap_name),
      `# signatures`     = sum(df$n_sigs, na.rm = TRUE),
      `Mean TAS`         = mean(df$tas, na.rm = TRUE),
      `Median TAS`       = median(df$tas, na.rm = TRUE),
      `TAS SD`           = sd(df$tas, na.rm = TRUE),
      `Cell lines`       = paste(sort(unique(df$cell_iname)), collapse = ", "),
      `Perturbation times` = paste(sort(unique(df$pert_time)), collapse = ", "),
      `Perturbation doses` = paste(sort(unique(df$pert_dose)), collapse = ", ")
    )
  })
  
  output$community_top_moas <- renderTable({
    df <- community_data() %>%
      dplyr::filter(!is.na(moa), moa != "")
    
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      dplyr::count(moa, name = "n_drugs", sort = TRUE) %>%
      dplyr::mutate(
        frac = n_drugs / sum(n_drugs)
      ) %>%
      dplyr::slice_head(n = 10)
  })
  
  output$community_top_atc <- renderTable({
    df <- community_data()
    
    # which ATC columns exist for this dataset?
    atc_cols <- intersect(c("level1", "level2", "level3", "level4"),
                          colnames(df))
    if (length(atc_cols) == 0) return(NULL)
    
    # long format: one row per (drug, level, atc_code)
    long_df <- df %>%
      dplyr::select(dplyr::all_of(atc_cols)) %>%
      tidyr::pivot_longer(
        cols      = everything(),
        names_to  = "level",
        values_to = "atc"
      ) %>%
      dplyr::filter(!is.na(atc), atc != "")
    
    if (nrow(long_df) == 0) return(NULL)
    
    # counts and fractions per level
    long_df %>%
      dplyr::count(level, atc, name = "n_drugs", sort = TRUE) %>%
      dplyr::group_by(level) %>%
      dplyr::mutate(
        frac = n_drugs / sum(n_drugs)
      ) %>%
      dplyr::slice_head(n = 10) %>%   # top 10 per level (tweak if you like)
      dplyr::ungroup() %>%
      dplyr::arrange(level, dplyr::desc(n_drugs))
  })
  
  
  output$community_drugs <- renderTable({
    df <- community_data()
    
    df %>%
      select(cmap_name, tas, pert_time, pert_dose) %>%
      arrange(desc(tas))
  })
  
  output$summary_title <- renderText({
    paste("Community", selected_community(), "Summary")
  })
  
  output$community_top_atc <- renderTable({
  df <- community_data()

  # which ATC columns exist for this dataset?
  atc_cols <- intersect(c("level1", "level2", "level3", "level4"),
                        colnames(df))
  if (length(atc_cols) == 0) return(NULL)

  # long format: one row per (drug, level, atc_code)
  long_df <- df %>%
    dplyr::select(dplyr::all_of(atc_cols)) %>%
    tidyr::pivot_longer(
      cols      = everything(),
      names_to  = "level",
      values_to = "atc"
    ) %>%
    dplyr::filter(!is.na(atc), atc != "")

  if (nrow(long_df) == 0) return(NULL)

  # counts and fractions per level
  long_df %>%
    dplyr::count(level, atc, name = "n_drugs", sort = TRUE) %>%
    dplyr::group_by(level) %>%
    dplyr::mutate(
      frac = n_drugs / sum(n_drugs)
    ) %>%
    dplyr::slice_head(n = 10) %>%   # top 10 per level (tweak if you like)
    dplyr::ungroup() %>%
    dplyr::arrange(level, dplyr::desc(n_drugs))
})

  output$edge_ridge_all <- renderPlot({
    # keep only edges where both endpoints are in the same community
    df <- edge_comm %>%
      dplyr::filter(!is.na(comm_src), !is.na(comm_dst), comm_src == comm_dst) %>%
      dplyr::transmute(
        community = comm_src,
        weight    = weight
      )
    
    if (nrow(df) == 0) return(NULL)
    
    # per-community summary: median + IQR + n_edges
    stats <- df %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        median_w = median(weight, na.rm = TRUE),
        q25      = quantile(weight, 0.25, na.rm = TRUE),
        q75      = quantile(weight, 0.75, na.rm = TRUE),
        n_edges  = dplyr::n(),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        # label shows community + number of edges
        community_label = paste0(community, " (n=", n_edges, ")")
      )
    
    # join stats back to each edge & order communities by median weight
    df2 <- df %>%
      dplyr::left_join(stats %>% dplyr::select(community, median_w, n_edges, community_label),
                       by = "community") %>%
      dplyr::mutate(
        community_label = stats::reorder(community_label, median_w)
      )
    
    # same ordering + labels for stats table (for vlines)
    stats2 <- stats %>%
      dplyr::mutate(
        community_label = stats::reorder(community_label, median_w)
      )
    
    ggplot(df2, aes(x = weight, y = community_label, fill = median_w)) +
      # main ridges
      ggridges::geom_density_ridges(
        alpha          = 0.8,
        scale          = 1.2,
        rel_min_height = 0.01,
        color          = "black",
        size           = 0.25
      ) +
      # # optional rug: comment out if it’s too dense
      # geom_rug(alpha = 0.05, sides = "b") +
      # # median (dashed) and IQR (dotted) lines
      # geom_vline(
      #   data = stats2,
      #   aes(xintercept = median_w),
      #   linetype = "dashed",
      #   alpha    = 0.35
      # ) +
      # geom_vline(
      #   data = stats2,
      #   aes(xintercept = q25),
      #   linetype = "dotted",
      #   alpha    = 0.25
      # ) +
      # geom_vline(
      #   data = stats2,
      #   aes(xintercept = q75),
      #   linetype = "dotted",
      #   alpha    = 0.25
      # ) +
      scale_fill_viridis_c(
        option = "A",
        name   = "Median\nedge weight"
      ) +
      labs(
        title = "Within-community signature-level edge weights (all communities)",
        x     = "Edge weight",
        y     = "Community (ordered by median weight)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text.y     = element_text(size = 7)
      )
  })
  
  output$edge_ridge_between <- renderPlot({
    # edges where endpoints are in *different* communities
    df <- edge_comm %>%
      dplyr::filter(!is.na(comm_src), !is.na(comm_dst), comm_src != comm_dst)
    
    if (nrow(df) == 0) return(NULL)
    
    # long format: each edge contributes once to each incident community
    df_long <- dplyr::bind_rows(
      df %>%
        dplyr::transmute(community = comm_src, partner = comm_dst, weight),
      df %>%
        dplyr::transmute(community = comm_dst, partner = comm_src, weight)
    )
    
    # per-community summary: median between weight + n_edges
    stats <- df_long %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        median_w = median(weight, na.rm = TRUE),
        n_edges  = dplyr::n(),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        community_label = paste0(community, " (n=", n_edges, ")")
      )
    
    # join stats back & order communities by median between weight
    df2 <- df_long %>%
      dplyr::left_join(
        stats %>% dplyr::select(community, median_w, community_label),
        by = "community"
      ) %>%
      dplyr::mutate(
        community_label = stats::reorder(community_label, median_w)
      )
    
    ggplot(df2, aes(x = weight, y = community_label, fill = median_w)) +
      ggridges::geom_density_ridges(
        alpha          = 0.8,
        scale          = 1.2,
        rel_min_height = 0.01,
        color          = "black",
        size           = 0.25
      ) +
      scale_fill_viridis_c(
        option = "C",
        name   = "Median\nedge weight"
      ) +
      labs(
        title    = "Between-community signature-level edge weights (all communities)",
        subtitle = "Edges connecting each community to other communities",
        x        = "Edge weight",
        y        = "Community (ordered by median between-community weight)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text.y     = element_text(size = 7)
      )
  })
  
  output$robustness_scatter <- renderPlot({
    # 1. Within-community stats
    within_stats <- edge_comm %>%
      dplyr::filter(!is.na(comm_src), !is.na(comm_dst), comm_src == comm_dst) %>%
      dplyr::group_by(community = comm_src) %>%
      dplyr::summarise(
        med_within = median(weight, na.rm = TRUE),
        n_within   = dplyr::n(),
        .groups    = "drop"
      )
    
    # 2. Between-community stats
    between_long <- edge_comm %>%
      dplyr::filter(!is.na(comm_src), !is.na(comm_dst), comm_src != comm_dst) %>%
      dplyr::bind_rows(
        # each edge contributes once for each endpoint community
        dplyr::transmute(., community = comm_src, weight),
        dplyr::transmute(., community = comm_dst, weight)
      )
    
    between_stats <- between_long %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        med_between = median(weight, na.rm = TRUE),
        n_between   = dplyr::n(),
        .groups     = "drop"
      )
    
    # 3. Combine + robustness score
    stats <- dplyr::full_join(within_stats, between_stats, by = "community") %>%
      dplyr::mutate(
        robustness = med_within - med_between,
        community  = as.character(community)
      )
    
    # drop communities missing either side
    stats <- stats %>%
      dplyr::filter(!is.na(med_within), !is.na(med_between))
    
    if (nrow(stats) == 0) return(NULL)
    
    ggplot(stats, aes(x = med_within, y = med_between)) +
      # diagonal: where within == between
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
      geom_point(
        aes(
          color = robustness,
          size  = n_within
        ),
        alpha = 0.9
      ) +
      ggrepel::geom_text_repel(
        aes(label = community),
        size        = 3,
        max.overlaps = 50,
        alpha       = 0.8
      ) +
      scale_color_viridis_c(
        option = "C",
        name   = "Robustness\n(within - between)"
      ) +
      scale_size_continuous(
        name = "# within edges",
        range = c(2, 8)
      ) +
      labs(
        title    = "Community robustness: within vs between edge weights",
        subtitle = "Good communities lie below the diagonal (within > between)",
        x        = "Median within-community edge weight",
        y        = "Median between-community edge weight"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right"
      )
  })
  
  # --- Compute per-community within and between statistics --- #
  comm_stats <- reactive({
    
    req(edge_comm)
    
    # 1) WITHIN-community edge stats
    within_stats <- edge_comm %>%
      dplyr::filter(
        !is.na(comm_src),
        !is.na(comm_dst),
        comm_src == comm_dst  # within
      ) %>%
      dplyr::group_by(community = comm_src) %>%
      dplyr::summarise(
        med_within = median(weight, na.rm = TRUE),
        mean_within = mean(weight, na.rm = TRUE),
        sd_within   = sd(weight, na.rm = TRUE),
        n_within    = n(),
        .groups = "drop"
      )
    
    # 2) BETWEEN-community edges
    # assign each endpoint community separately
    between_long <- edge_comm %>%
      dplyr::filter(
        !is.na(comm_src),
        !is.na(comm_dst),
        comm_src != comm_dst
      ) %>%
      dplyr::bind_rows(
        dplyr::transmute(., community = comm_src, weight),
        dplyr::transmute(., community = comm_dst, weight)
      )
    
    between_stats <- between_long %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        med_between = median(weight, na.rm = TRUE),
        mean_between = mean(weight, na.rm = TRUE),
        sd_between   = sd(weight, na.rm = TRUE),
        n_between    = n(),
        .groups = "drop"
      )
    
    # 3) Combine stats
    stats <- dplyr::full_join(within_stats, between_stats, by = "community") %>%
      dplyr::mutate(
        robustness = med_within - med_between,
        community  = as.character(community)
      ) %>%
      dplyr::filter(!is.na(med_within), !is.na(med_between))
    
    stats
  })
  
  output$robustness_dumbbell <- renderPlot({
    
    stats <- comm_stats()
    req(nrow(stats) > 0)
    
    # Order by robustness
    stats <- stats %>% arrange(desc(robustness))
    
    ggplot(stats, aes(y = reorder(community, robustness))) +
      
      # Dumbbell line
      geom_segment(
        aes(
          x    = med_between,
          xend = med_within,
          yend = community
        ),
        color = "grey70",
        size  = 1.1,
        alpha = 0.6
      ) +
      
      # Left point (between) — same color scale but faded
      geom_point(
        aes(x = med_between, color = robustness),
        size = 3,
        alpha = 0.55
      ) +
      
      # Right point (within)
      geom_point(
        aes(x = med_within, color = robustness, size = log10(n_within + 1)),
        alpha = 0.95
      ) +
      
      scale_color_viridis_c(
        option = "C",
        name   = "Robustness\n(within - between)"
      ) +
      
      scale_size_continuous(
        name = "# within edges\n(log10)",
        range = c(1, 8)
      ) +
      
      labs(
        title = "Community robustness dumbbell: within vs between edge weights",
        subtitle = "Left = between-community median weight; Right = within-community median weight",
        x = "Median edge weight",
        y = "Community (ordered by robustness)"
      ) +
      
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        panel.grid.major.y = element_blank()
      ) + 
      
      ggrepel::geom_text_repel(
        aes(
          x = med_within,
          label = community
        ),
        size = 3,
        nudge_x = 0.02,
        direction = "y",
        segment.alpha = 0.1,
        color = "black"
      )
    
  })
  
  # ---- MoA ↔ ATC cross-ontology data (drug-level) ----
  moa_atc_data <- reactive({
    df <- filtered_nodes()
    
    # clean MoA: trim whitespace & drop unlabeled
    df <- df %>%
      dplyr::mutate(
        moa_clean = trimws(moa)
      ) %>%
      dplyr::filter(
        !is.na(moa_clean),
        moa_clean != "\"\"",
        moa_clean != "NA"        # just in case string "NA" appears
      )
    
    # stop if nothing left
    if (nrow(df) == 0) return(NULL)
    
    lvl <- input$atc_level  # "level1"..."level4"
    if (!lvl %in% colnames(df)) return(NULL)
    
    df <- df %>%
      dplyr::filter(!is.na(.data[[lvl]]), .data[[lvl]] != "")
    
    if (nrow(df) == 0) return(NULL)
    
    # distinct drug–MoA–ATC combos, then count drugs
    flows <- df %>%
      dplyr::distinct(cmap_name, moa = moa_clean, atc = .data[[lvl]]) %>%
      dplyr::count(moa, atc, name = "n_drugs")
    
    # keep only reasonably supported flows
    min_drugs <- 3
    flows <- flows %>% dplyr::filter(n_drugs >= min_drugs)
    if (nrow(flows) == 0) return(NULL)
    
    # limit to top MoAs and top ATC classes so the plot stays readable
    top_moas <- flows %>%
      dplyr::group_by(moa) %>%
      dplyr::summarise(total = sum(n_drugs), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::pull(moa)
    
    top_atc <- flows %>%
      dplyr::group_by(atc) %>%
      dplyr::summarise(total = sum(n_drugs), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::pull(atc)
    
    flows %>%
      dplyr::filter(moa %in% top_moas, atc %in% top_atc)
  })
  
  output$moa_atc_alluvial <- renderPlot({
    flows <- moa_atc_data()
    req(!is.null(flows), nrow(flows) > 0)
    
    # human-readable label for the ATC level
    atc_label <- switch(
      input$atc_level,
      "level1" = "ATC Level 1",
      "level2" = "ATC Level 2",
      "level3" = "ATC Level 3",
      "level4" = "ATC Level 4",
      input$atc_level
    )
    
    ggplot(
      flows,
      aes(
        axis1 = moa,
        axis2 = atc,
        y     = n_drugs
      )
    ) +
      ggalluvial::geom_alluvium(
        aes(fill = moa),
        width = 1/12,
        alpha = 0.8,
        knot.pos = 0.4
      ) +
      ggalluvial::geom_stratum(
        width = 1/6,
        fill  = "grey95",
        color = "grey50"
      ) +
      geom_text(
        stat  = "stratum",
        aes(label = after_stat(stratum)),
        size  = 3
      ) +
      scale_x_discrete(
        limits = c("MoA", atc_label),
        expand = c(.1, .05)
      ) +
      labs(
        title    = "MoA \u2194 ATC cross-ontology mapping",
        subtitle = paste0(
          "Flows weighted by # of distinct drugs; showing top MoAs and ",
          atc_label, " classes after current filters"
        ),
        x = NULL,
        y = "# Drugs"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position    = "none"
      )
  })
  
  # --- Per-community MoA purity (after filters) --- #
  moa_purity_data <- reactive({
    df <- filtered_nodes() %>%
      dplyr::filter(!is.na(moa), moa != "\"\"")
    
    if (nrow(df) == 0) return(tibble::tibble())
    
    # count MoAs within each community
    by_moa <- df %>%
      dplyr::count(community, moa, name = "n_moa")
    
    # total drugs with MoA per community
    totals <- by_moa %>%
      dplyr::group_by(community) %>%
      dplyr::summarise(
        n_total = sum(n_moa),
        .groups = "drop"
      )
    
    # pick the dominant MoA per community and compute purity
    purity <- by_moa %>%
      dplyr::group_by(community) %>%
      dplyr::slice_max(n_moa, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(totals, by = "community") %>%
      dplyr::mutate(
        purity   = n_moa / n_total,
        community = as.character(community)
      )
    
    purity
  })
  
  output$moa_entropy <- renderPlot({
    stats <- moa_entropy_data()
    if (is.null(stats) || nrow(stats) == 0) return(NULL)
    
    # order communities from most pure (low entropy) to most mixed (high)
    stats <- stats %>%
      dplyr::arrange(entropy) %>%
      dplyr::mutate(
        community_factor = reorder(as.character(community), entropy),
        label = sprintf(
          "%s (%.0f%%, n=%d)",
          dominant_moa,
          100 * dom_frac,
          n_drugs
        )
      )
    
    ggplot(
      stats,
      aes(
        x     = entropy,
        y     = community_factor,
        fill  = entropy
      )
    ) +
      geom_col(alpha = 0.9) +
      geom_text(
        aes(label = label),
        hjust = -0.05,
        size  = 3
      ) +
      scale_x_continuous(
        limits = c(0, 1.05),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::number_format(accuracy = 0.1)
      ) +
      scale_fill_viridis_c(
        option = "C",
        name   = "MoA entropy\n(0 = pure, 1 = mixed)"
      ) +
      labs(
        title    = "Community MoA entropy",
        subtitle = "Entropy computed over *labeled* MoAs only (unlabeled drugs are ignored)",
        x        = "MoA entropy",
        y        = "Community"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position      = "right",
        panel.grid.major.y   = element_blank(),
        axis.text.y          = element_text(size = 7)
      )
  })
  
  
  
  
  
}