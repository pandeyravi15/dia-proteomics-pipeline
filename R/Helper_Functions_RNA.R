# ==============================================================================
# Helper_Functions.R
# Brain RNA-seq Analysis — Core Utility Functions
#
# Description:
#   Shared functions for transcriptomics QC, DESeq2 differential expression,
#   KEGG enrichment, co-expression module correlation, AD subtype correlation,
#   and BioDomain ORA/GSEA. Designed for any case-control RNA-seq study with
#   sex and age stratification.
#
# Usage:
#   source("R/Helper_Functions.R")
#   All user-configurable settings (paths, organism, thresholds) should be
#   defined in the USER CONFIGURATION block in the main RMD notebook.
#
# Changelog:
#   - Removed all hardcoded Box paths — all reference data loaded via notebook config
#   - organism_db is now a parameter, not hardcoded to org.Mm.eg.db
#   - kegg_enrichment_analysis() organism is now a parameter (default: "mmu")
#   - enrichORA_BD() and enrichGSEA_BD() accept org_db and bio_ref as arguments
#   - run_lm_subset() generalizes LM across any formula and age subset
#   - All functions have NULL guards for empty results
#
# Author: Ravi S Pandey
# ==============================================================================


# ==============================================================================
# 1. SETUP & DEPENDENCIES
# ==============================================================================

needed.packages <- c(
  "ComplexHeatmap", "AnnotationDbi", "circlize", "corrplot", "cowplot",
  "clusterProfiler", "dplyr", "DESeq2", "DOSE", "data.table",
  "EnhancedVolcano", "enrichplot", "forcats", "ggplot2", "gridExtra",
  "ggpubr", "ggnewscale", "ggplotify", "ggrepel", "ggh4x", "Hmisc",
  "janitor", "lubridate", "limma", "purrr", "PCAtools",
  "ReactomePA", "RColorBrewer", "stringr", "scales", "sva",
  "tidyr", "tibble", "UpSetR", "VennDiagram", "gt",
  "DEGreport", "gprofiler2", "readxl", "readr", "broom", "matrixStats"
)

suppressPackageStartupMessages({
  lapply(needed.packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0(
        "Required package '", pkg, "' is not installed.\n",
        "Install with: install.packages('", pkg, "') ",
        "or BiocManager::install('", pkg, "')"
      ))
    }
    library(pkg, character.only = TRUE)
  })
})


# ==============================================================================
# 2. UTILITY: PCA
# ==============================================================================

#' Generate PCA Data Frame from Expression Matrix
#'
#' @param object   Matrix or data frame of normalized counts/TPM (genes as rows)
#' @param metadata Data frame with sample information
#' @param intgroup Character vector of columns in metadata to include
#' @param ntop     Number of top variable genes to use for PCA
#' @param returnData Logical; if TRUE returns data frame with percentVar attribute
#' @return Data frame of PC scores with metadata columns attached
plotPCA.df <- function(object, metadata, intgroup = "condition",
                       ntop = 500, returnData = TRUE) {

  rv     <- matrixStats::rowVars(as.matrix(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca    <- prcomp(t(object[select, ]))

  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (!all(intgroup %in% names(metadata)))
    stop("'intgroup' must specify valid columns in the metadata.")

  n_pcs <- min(10, ncol(pca$x))
  d     <- as.data.frame(pca$x[, 1:n_pcs])
  d     <- cbind(d, metadata[, intgroup, drop = FALSE])
  d$name <- rownames(metadata)

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:n_pcs]
    return(d)
  }
}


# ==============================================================================
# 3. UTILITY: GENE ID MAPPING
# ==============================================================================

#' Map Gene IDs Between Types Using an OrgDb Object
#'
#' Strips Ensembl version suffixes before mapping.
#' Works with any organism — pass org_db as argument rather than hardcoding.
#'
#' @param x           Vector of gene IDs, or a data frame/matrix with rownames as IDs
#' @param input_type  Input key type (e.g., "ENSEMBL", "ENTREZID")
#' @param output_type Output key type (e.g., "SYMBOL", "ENTREZID")
#' @param org_db      OrgDb object (e.g., org.Mm.eg.db). Defaults to org.Mm.eg.db.
#' @return Character vector of mapped IDs (NA where mapping fails)
map_gene_ids <- function(x, input_type = "ENSEMBL", output_type = "SYMBOL",
                         org_db = org.Mm.eg.db) {

  gene_keys <- if (is.vector(x)) x else rownames(x)
  if (is.null(gene_keys)) stop("Input must have rownames or be a character vector.")

  # Strip Ensembl version suffixes (e.g., ENSMUSG000...1.2 → ENSMUSG000...1)
  gene_keys <- gsub("\\..*$", "", as.character(gene_keys))

  mapped_ids <- tryCatch(
    AnnotationDbi::mapIds(
      org_db,
      keys      = gene_keys,
      column    = output_type,
      keytype   = input_type,
      multiVals = "first"
    ),
    error = function(e) rep(NA_character_, length(gene_keys))
  )

  na_count <- sum(is.na(mapped_ids))
  message(sprintf("Mapped %d/%d IDs successfully (%d NA).",
                  length(gene_keys) - na_count, length(gene_keys), na_count))

  return(as.character(mapped_ids))
}


# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION: DESeq2
# ==============================================================================

#' Run DESeq2 for a Single Pairwise Comparison
#'
#' Subsets counts and metadata to two groups, builds a DESeqDataSet, runs
#' DESeq2, and annotates results with gene symbols and Entrez IDs.
#'
#' @param raw_counts  Matrix/data frame of raw integer counts (genes × samples)
#' @param metadata    Data frame where rownames match colnames of raw_counts
#' @param design_col  Column name in metadata to use as the design factor
#' @param min_counts  Minimum total row sum to retain a gene (default: 10)
#' @param alpha       Significance threshold for results() (default: 0.05)
#' @param org_db      OrgDb object for annotation (default: org.Mm.eg.db)
#' @return Data frame of DESeq2 results with symbol and EntrezGene columns
analyze_deg_sample_type <- function(raw_counts, metadata, design_col,
                                    min_counts = 10, alpha = 0.05,
                                    org_db = org.Mm.eg.db) {

  common_samples <- intersect(colnames(raw_counts), rownames(metadata))
  if (length(common_samples) == 0)
    stop("No matching sample IDs between raw_counts columns and metadata rows.")

  design_formula <- as.formula(paste0("~", design_col))
  counts_sub     <- as.matrix(raw_counts[, common_samples])
  meta_sub       <- metadata[common_samples, , drop = FALSE]

  dds  <- DESeq2::DESeqDataSetFromMatrix(counts_sub, meta_sub, design_formula)
  dds  <- dds[rowSums(DESeq2::counts(dds)) >= min_counts, ]
  dds  <- DESeq2::DESeq(dds, parallel = TRUE)
  res  <- DESeq2::results(dds, alpha = alpha)

  res$symbol     <- map_gene_ids(res, "ENSEMBL", "SYMBOL",    org_db)
  res$EntrezGene <- map_gene_ids(res, "ENSEMBL", "ENTREZID",  org_db)

  final_res <- as.data.frame(res)
  col_order  <- c("symbol", "EntrezGene",
                  setdiff(names(final_res), c("symbol", "EntrezGene")))
  return(final_res[, col_order])
}


#' Run DESeq2 Across All Pairwise Comparisons
#'
#' Iterates over a comparison table, subsets data, and calls
#' analyze_deg_sample_type() for each pair. Returns combined list and data frame.
#'
#' @param counts          Raw count matrix
#' @param meta            Metadata with a Group column
#' @param comp_table      Data frame with columns: control, case
#' @param target_col      Column name used as DESeq2 design factor
#' @param experiment_name String label for this experiment
#' @param org_db          OrgDb object for annotation
#' @return List with $df (combined data frame) and $list (named list per comparison)
run_experiment_analysis <- function(counts, meta, comp_table,
                                    target_col, experiment_name,
                                    org_db = org.Mm.eg.db) {

  if (nrow(comp_table) == 0)
    stop("comp_table is empty — no valid comparisons to run.")

  exp_list <- list()
  exp_df   <- data.frame()

  for (i in 1:nrow(comp_table)) {
    ctrl_grp <- comp_table$control[i]
    case_grp <- comp_table$case[i]

    message("Running DESeq2: ", case_grp, " vs ", ctrl_grp)

    meta_sub             <- meta %>% dplyr::filter(Group %in% c(ctrl_grp, case_grp))
    rownames(meta_sub)   <- meta_sub$Names

    # Skip if either group has fewer than 3 samples
    group_counts <- table(as.character(meta_sub$Group))
    if (any(group_counts < 3)) {
      message("  Skipping — insufficient samples (n < 3) in at least one group.")
      next
    }

    res <- tryCatch(
      analyze_deg_sample_type(counts, meta_sub, design_col = target_col,
                               org_db = org_db),
      error = function(e) {
        message("  DESeq2 failed: ", e$message)
        return(NULL)
      }
    )
    if (is.null(res)) next

    # Build readable label from case group name
    case_parts <- unlist(strsplit(case_grp, "_"))
    label      <- paste0(case_parts[1], "(", case_parts[2], "_", case_parts[3], ")")

    exp_list[[label]] <- res
    exp_df <- rbind(exp_df,
                    res %>% dplyr::mutate(group = label,
                                          ctrl_group  = ctrl_grp,
                                          Experiment  = experiment_name))
  }

  message("DESeq2 complete: ", length(exp_list), " of ", nrow(comp_table),
          " comparisons succeeded.")
  return(list(df = exp_df, list = exp_list))
}


# ==============================================================================
# 5. FUNCTIONAL ENRICHMENT: KEGG
# ==============================================================================

#' Run KEGG Enrichment on Multiple Groups
#'
#' @param groups       Vector of group names to analyze
#' @param deg_results  Data frame with columns: group, log2FoldChange, EntrezGene
#' @param top_n        Number of categories for dotplot (default: 10)
#' @param organism     KEGG organism code (default: "mmu"; "hsa" for human)
#' @param pval_cutoff  p-value cutoff for enrichment (default: 0.05)
kegg_enrichment_analysis <- function(groups, deg_results, top_n = 10,
                                      organism = "mmu", pval_cutoff = 0.05) {

  get_gene_list <- function(grp_name, direction = "up") {
    filter_expr <- if (direction == "up") deg_results$log2FoldChange > 0 else
      deg_results$log2FoldChange < 0
    deg_results %>%
      dplyr::filter(group == grp_name & filter_expr) %>%
      dplyr::pull(EntrezGene) %>%
      na.omit() %>% unique() %>% as.character()
  }

  run_direction <- function(dir_label, dir) {
    gene_bundles <- lapply(groups, get_gene_list, direction = dir)
    names(gene_bundles) <- groups

    if (sum(lengths(gene_bundles)) == 0) {
      message("No genes found for direction: ", dir_label)
      return(invisible(NULL))
    }

    # Drop groups with too few genes
    gene_bundles <- gene_bundles[sapply(gene_bundles, length) >= 3]
    if (length(gene_bundles) == 0) {
      message("No groups with >= 3 genes for direction: ", dir_label)
      return(invisible(NULL))
    }

    enr_res <- tryCatch(
      clusterProfiler::compareCluster(gene_bundles, fun = "enrichKEGG",
                                       pvalueCutoff = pval_cutoff,
                                       organism = organism),
      error = function(e) {
        message("KEGG error (", dir_label, "): ", e$message)
        return(NULL)
      }
    )

    if (is.null(enr_res) || nrow(as.data.frame(enr_res)) == 0) {
      message("No significant KEGG pathways for ", dir_label)
      return(invisible(NULL))
    }

    # Clean species suffix from pathway descriptions
    enr_res@compareClusterResult$Description <- gsub(
      " - .*$", "", enr_res@compareClusterResult$Description
    )

    print(
      clusterProfiler::dotplot(enr_res, showCategory = top_n,
                                label_format = 60) +
        ggplot2::ggtitle(paste("KEGG Pathway Enrichment:", dir_label)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
          plot.title  = ggplot2::element_text(face = "bold", size = 14)
        ) +
        ggplot2::xlab("Comparison Groups")
    )
  }

  run_direction("Upregulated",   "up")
  run_direction("Downregulated", "down")
}


# ==============================================================================
# 6. CORRELATION: MOUSE vs. HUMAN AMP-AD MODULES (logFC)
# ==============================================================================

#' Correlate Mouse Fold-Changes with Human AMP-AD Module Data
#'
#' @param mouse_res   Data frame with columns: symbol (gene), log2FoldChange, group
#' @param human_data  Data frame with columns: Gene, ampad_fc, module
#' @param module_meta Data frame mapping module to cluster (cluster, cluster_label)
#' @param mod_order   Factor levels for modules (x-axis order)
#' @param sex_order   Factor levels for groups (y-axis order)
#' @return Data frame ready for magora_corrplot_colored()
correlation_analysis <- function(mouse_res, human_data, module_meta,
                                  mod_order, sex_order) {

  combined <- mouse_res %>%
    dplyr::rename(Gene = symbol) %>%
    dplyr::inner_join(human_data, by = "Gene", relationship = "many-to-many") %>%
    dplyr::select(module, group, Gene, log2FoldChange, ampad_fc)

  cor_results <- combined %>%
    dplyr::group_by(module, group) %>%
    tidyr::nest(data = c(Gene, log2FoldChange, ampad_fc)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, purrr::possibly(
        ~cor.test(.x$log2FoldChange, .x$ampad_fc, method = "pearson"),
        otherwise = NULL
      )),
      correlation = purrr::map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value     = purrr::map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    dplyr::filter(!is.na(correlation)) %>%
    dplyr::ungroup()

  cor_results %>%
    dplyr::mutate(significant = p_value < 0.05) %>%
    dplyr::left_join(module_meta, by = "module") %>%
    dplyr::mutate(
      module    = factor(module, levels = mod_order),
      model_sex = factor(group,  levels = sex_order),
      model_sex = forcats::fct_rev(model_sex)
    ) %>%
    dplyr::select(cluster, cluster_label, module, model_sex,
                  correlation, p_value, significant) %>%
    dplyr::arrange(cluster, model_sex)
}


# ==============================================================================
# 7. VISUALIZATION: MODULE CORRELATION (logFC)
# ==============================================================================

#' AMP-AD Module Correlation Dot Plot with Colored Facet Strips
#'
#' @param data Output from correlation_analysis()
magora_corrplot_colored <- function(data) {

  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100

  p <- ggplot2::ggplot(data, ggplot2::aes(x = module, y = model_sex)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation)), alpha = 0.8) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      ggplot2::aes(size = abs(correlation)),
      color = "black", shape = 0, size = 10, stroke = 1.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(range = c(1, 8), limits = c(0, ran),
                                   guide = "none") +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", mid = "white", high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE, frame.colour = "black")
    ) +
    ggplot2::facet_grid(cols = dplyr::vars(cluster_label),
                        scales = "free_x", space = "free_x", switch = "y") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.x      = ggplot2::element_text(size = 14, face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "grey95", colour = NA),
      axis.text.x       = ggplot2::element_text(angle = 90, hjust = 0, size = 14),
      axis.text.y       = ggplot2::element_text(size = 15),
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "right"
    )

  # Color facet strip backgrounds
  fills  <- c("darkorange3", "chartreuse3", "deepskyblue2", "turquoise", "deeppink2")
  g      <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  stripr <- which(grepl("strip-t", g$layout$name))
  if (length(stripr) > length(fills))
    fills <- rep(fills, length.out = length(stripr))

  for (i in seq_along(stripr)) {
    j <- which(grepl("rect", g$grobs[[stripr[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[stripr[i]]]$grobs[[1]]$children[[j]]$gp$fill <- fills[i]
  }

  cowplot::ggdraw(g)
}


# ==============================================================================
# 8. CORRELATION: MOUSE vs. HUMAN AMP-AD MODULES (LM Estimates)
# ==============================================================================

#' Correlate LM Estimates with Human AMP-AD Module Data
#'
#' @param data        Data frame with columns: Gene, Variant, value (LM estimates)
#' @param human_data  Data frame with columns: Gene, ampad_fc, module
#' @param module_meta Mapping of module to cluster
#' @param mod_order   Factor levels for modules
#' @param var_order   Factor levels for variants
#' @return Data frame ready for variant_corrplot_ggh4x()
corr_function_lm <- function(data, human_data, module_meta, mod_order, var_order) {

  lm_vs_human <- data %>%
    dplyr::inner_join(human_data, by = "Gene") %>%
    dplyr::group_by(module, Variant) %>%
    dplyr::filter(dplyr::n() > 2) %>%
    tidyr::nest(data_subset = c(Gene, value, ampad_fc)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data_subset, purrr::possibly(
        ~cor.test(.x$value, .x$ampad_fc, method = "pearson"), otherwise = NULL
      )),
      correlation = purrr::map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value     = purrr::map_dbl(cor_test, ~if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    dplyr::filter(!is.na(correlation)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-cor_test, -data_subset)

  lm_vs_human %>%
    dplyr::mutate(
      significant = p_value < 0.05,
      age_group   = "All Months"
    ) %>%
    dplyr::left_join(module_meta, by = "module") %>%
    dplyr::mutate(
      Variant = factor(Variant, levels = var_order),
      Variant = forcats::fct_rev(Variant),
      module  = factor(module, levels = mod_order)
    ) %>%
    dplyr::select(cluster, cluster_label, module, Variant,
                  age_group, correlation, p_value, significant) %>%
    dplyr::arrange(cluster, Variant)
}


# ==============================================================================
# 9. VISUALIZATION: LM MODULE CORRELATION
# ==============================================================================

#' LM Variant × Module Correlation Dot Plot (ggh4x faceting)
#'
#' Requires a `Background` column in data for row faceting.
#' @param data Data frame from corr_function_lm() with Background column added
#' @param ran  Numeric. Color scale limit
variant_corrplot_ggh4x <- function(data, ran) {

  ggplot2::ggplot(data, ggplot2::aes(x = module, y = Variant)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation))) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      shape = 0, size = 9, colour = "black", stroke = 1
    ) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(Background),
      cols = ggplot2::vars(cluster_label),
      scales = "free", space = "free", switch = "y",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = c("darkorange3", "chartreuse3", "deepskyblue2",
                   "turquoise", "deeppink2")
        )
      )
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      strip.text.x       = ggplot2::element_text(size = 11, face = "bold",
                                                   color = "black"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 14),
      strip.background.y = ggplot2::element_rect(fill = "grey95"),
      axis.text.x        = ggplot2::element_text(angle = 90, hjust = 0, size = 14),
      axis.text.y        = ggplot2::element_text(angle = 0, hjust = 0, size = 14),
      panel.grid         = ggplot2::element_blank(),
      legend.position    = "right",
      panel.spacing      = ggplot2::unit(0.2, "lines")
    )
}


# ==============================================================================
# 10. CORRELATION: MOUSE vs. HUMAN AD SUBTYPES
# ==============================================================================

#' Correlate DESeq2 logFC with Human AD Subtype Signatures
#'
#' Handles both Nikhil et al. (with cohort metadata) and Neff et al.
#' (without cohort metadata, using a custom label).
#'
#' @param load2_fc       DESeq2 data frame with symbol, log2FoldChange, group
#' @param target_data    Reference data frame with Gene, ampad_fc, subtype
#' @param cohort_data    Optional cohort metadata (subtype → cluster_label)
#' @param custom_label   String label if cohort_data is NULL
#' @param subtypes_levels Factor order for subtypes
#' @param sex_levels      Factor order for groups (y-axis)
#' @return Data frame ready for subtype_FC_corrplot_colored()
run_subtype_correlation <- function(load2_fc, target_data,
                                     subtypes_levels = NULL,
                                     sex_levels      = NULL,
                                     cohort_data     = NULL,
                                     custom_label    = "MSBB_PHG_Neff") {

  results <- load2_fc %>%
    dplyr::rename(Gene = dplyr::any_of(c("symbol", "Gene"))) %>%
    dplyr::inner_join(target_data, by = "Gene") %>%
    dplyr::select(subtype, group, Gene, log2FoldChange, ampad_fc) %>%
    dplyr::group_by(subtype, group) %>%
    tidyr::nest(data = c(Gene, log2FoldChange, ampad_fc)) %>%
    dplyr::mutate(
      cor_test = purrr::map(data, purrr::possibly(
        ~stats::cor.test(.x$log2FoldChange, .x$ampad_fc, method = "pearson"),
        otherwise = NULL
      )),
      correlation = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value     = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$p.value),
      significant = p_value < 0.05
    ) %>%
    dplyr::filter(!is.na(correlation)) %>%
    dplyr::ungroup()

  if (!is.null(cohort_data)) {
    results <- results %>% dplyr::left_join(cohort_data, by = "subtype")
  } else {
    results <- results %>%
      dplyr::mutate(cluster = "Unknown", cluster_label = custom_label)
  }

  if (!is.null(sex_levels))
    results <- results %>%
      dplyr::mutate(
        model_sex = factor(group, levels = sex_levels),
        model_sex = forcats::fct_rev(model_sex)
      )

  if (!is.null(subtypes_levels))
    results$subtype <- factor(results$subtype, levels = subtypes_levels)

  return(results)
}


# ==============================================================================
# 11. VISUALIZATION: SUBTYPE CORRELATION
# ==============================================================================

#' Subtype Correlation Dot Plot with Colored Strip Backgrounds
#'
#' @param data Output from run_subtype_correlation()
#' @param ran  Numeric. Absolute maximum for correlation scale
subtype_FC_corrplot_colored <- function(data, ran = 1.0) {

  facet_colors    <- c("darkorange3", "chartreuse3", "deepskyblue2",
                        "turquoise", "deeppink2")
  n_labels        <- length(unique(data$cluster_label))
  strip_fills     <- rep(facet_colors, length.out = n_labels)

  ggplot2::ggplot(data, ggplot2::aes(x = subtype, y = model_sex)) +
    ggplot2::geom_tile(colour = "gray92", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation))) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      shape = 0, size = 10, stroke = 1, color = "black"
    ) +
    ggh4x::facet_grid2(
      cols  = ggplot2::vars(cluster_label),
      scales = "free", space = "free",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = strip_fills)
      )
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(guide = "none", range = c(1, 8),
                                   limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E",
      midpoint = 0, limits = c(-ran, ran), name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE, title.position = "top",
                                       title.hjust = 0.5)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      strip.text.x      = ggplot2::element_text(size = 11, face = "bold",
                                                  color = "black"),
      axis.text.x       = ggplot2::element_text(size = 10, angle = 90, hjust = 0),
      axis.text.y       = ggplot2::element_text(size = 10, face = "bold"),
      panel.spacing     = ggplot2::unit(0.3, "lines"),
      panel.border      = ggplot2::element_rect(color = "gray80", fill = NA),
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "bottom",
      legend.key.width  = ggplot2::unit(2, "cm")
    )
}


# ==============================================================================
# 12. BIODOMAIN CORRELATION: LM ESTIMATES
# ==============================================================================

#' Correlate LM Estimates with BioDomain Reference Fold-Changes
#'
#' @param data           Data frame with columns: Gene, Variant, value
#' @param ref            BioDomain reference with columns: Gene, Term, ampad_fc
#' @param sig            p-value threshold for significance flag
#' @param age_group_val  String label for this age group (used in output)
#' @param variant_levels Factor order for Variants
#' @param domain_levels  Factor order for BioDomains
#' @return Data frame ready for variant_corrplot_biodom()
corr_function_lm_biodomain <- function(data, ref, sig, age_group_val,
                                        variant_levels, domain_levels) {

  results <- data %>%
    dplyr::inner_join(ref, by = "Gene") %>%
    dplyr::filter(is.finite(value), is.finite(ampad_fc)) %>%
    dplyr::group_by(Term, Variant) %>%
    dplyr::filter(dplyr::n() >= 3, sd(value) > 0, sd(ampad_fc) > 0) %>%
    tidyr::nest(data = c(Gene, value, ampad_fc)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, purrr::possibly(
        ~stats::cor.test(.x$value, .x$ampad_fc, method = "pearson"),
        otherwise = NULL
      )),
      correlation = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value     = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    dplyr::filter(!is.na(correlation)) %>%
    dplyr::ungroup()

  results %>%
    dplyr::mutate(
      significant = p_value < sig,
      age_group   = age_group_val,
      Variant     = factor(Variant, levels = variant_levels),
      Variant     = forcats::fct_rev(Variant),
      module      = factor(Term, levels = domain_levels)
    ) %>%
    dplyr::select(biodom = Term, Variant, age_group,
                  correlation, p_value, significant, module) %>%
    dplyr::arrange(module)
}


#' BioDomain Variant Correlation Dot Plot
#'
#' Requires `group` and `Name` columns added by the notebook before calling.
#' @param data Data frame from corr_function_lm_biodomain() with group and Name columns
#' @param ran  Color scale limit
variant_corrplot_biodom <- function(data, ran) {

  ggplot2::ggplot(data, ggplot2::aes(x = biodom, y = Variant)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation)), alpha = 0.85) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      ggplot2::aes(size = abs(correlation)),
      color = "black", shape = 0, size = 10, stroke = 1.1,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(range = c(1, 8), limits = c(0, ran),
                                   guide = "none") +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", mid = "white", high = "#164B6E",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE, frame.colour = "black")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::facet_grid(
      rows = dplyr::vars(group), cols = dplyr::vars(Name),
      scales = "free", space = "free", switch = "y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.x       = ggplot2::element_text(size = 13, face = "bold"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      strip.background   = ggplot2::element_rect(fill = "grey96", colour = NA),
      axis.ticks         = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 90, hjust = 0, size = 12),
      axis.text.y        = ggplot2::element_text(size = 12),
      panel.grid         = ggplot2::element_blank(),
      legend.position    = "right"
    )
}


# ==============================================================================
# 13. SUB-BIODOMAIN CORRELATION: LM ESTIMATES
# ==============================================================================

#' Correlate LM Estimates with Sub-BioDomain Reference Fold-Changes
#'
#' @param data           Data frame with columns: Gene, Variant, value
#' @param ref            Sub-BioDomain reference with Gene, Term, ampad_fc, Biodomain
#' @param sig            p-value threshold
#' @param age_group_val  String label for this age group
#' @param variant_levels Factor order for Variants
#' @param domain_levels  Factor order for broad BioDomains
#' @return Data frame ready for variant_corrplot_subbiodom_portrait()
corr_function_lm_subbiodomain <- function(data, ref, sig, age_group_val,
                                           variant_levels, domain_levels) {

  results <- data %>%
    dplyr::rename(Gene = 1) %>%
    dplyr::inner_join(ref, by = "Gene") %>%
    dplyr::filter(is.finite(value), is.finite(ampad_fc)) %>%
    dplyr::group_by(Term, Variant) %>%
    dplyr::filter(dplyr::n() >= 3, sd(value) > 0, sd(ampad_fc) > 0) %>%
    tidyr::nest(data = c(Gene, value, ampad_fc)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, purrr::possibly(
        ~stats::cor.test(.x$value, .x$ampad_fc, method = "pearson"),
        otherwise = NULL
      )),
      correlation = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$estimate),
      p_value     = purrr::map_dbl(cor_test,
                                    ~if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    dplyr::filter(!is.na(correlation)) %>%
    dplyr::ungroup()

  results %>%
    dplyr::left_join(human_biodomains, by = "Biodomain") %>%
    dplyr::mutate(
      significant = p_value < sig,
      age_group   = age_group_val,
      Variant     = factor(Variant, levels = variant_levels),
      Variant     = forcats::fct_rev(Variant),
      module_fact = factor(Biodomain, levels = domain_levels)
    ) %>%
    dplyr::select(module = Biodomain, biodom = Term, Variant, age_group,
                  correlation, p_value, significant, Biodomain_label,
                  module_fact) %>%
    dplyr::arrange(module_fact)
}


#' Sub-BioDomain Portrait Correlation Dot Plot
#'
#' Requires `group` column added by the notebook before calling.
#' @param data Data frame from corr_function_lm_subbiodomain() with group column
#' @param ran  Color scale limit
variant_corrplot_subbiodom_portrait <- function(data, ran = 1.0) {

  ggplot2::ggplot(data, ggplot2::aes(y = biodom, x = Variant)) +
    ggplot2::geom_tile(colour = "gray92", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation))) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      shape = 0, size = 10, stroke = 0.8, color = "black"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(guide = "none", range = c(1, 8),
                                   limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E", midpoint = 0,
      limits = c(-ran, ran), name = "Correlation", na.value = "gray95",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Biodomain_label), cols = ggplot2::vars(group),
      scales = "free", space = "free", switch = "y"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.x       = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      strip.background   = ggplot2::element_rect(fill = "gray95", color = "gray80"),
      axis.text.x        = ggplot2::element_text(angle = 90, hjust = 0, size = 11),
      axis.text.y        = ggplot2::element_text(size = 10),
      panel.spacing      = ggplot2::unit(0.1, "lines"),
      panel.border       = ggplot2::element_rect(color = "gray80", fill = NA),
      panel.grid         = ggplot2::element_blank(),
      legend.position    = "right"
    )
}


#' Sub-BioDomain Landscape Correlation Dot Plot
#'
#' Same as portrait but with axes transposed (biodom on x, Variant on y).
#' @param data Filtered sub-domain correlation data frame with group column
#' @param ran  Color scale limit
variant_corrplot_subbiodom_land <- function(data, ran = 1.0) {

  ggplot2::ggplot(data, ggplot2::aes(x = biodom, y = Variant)) +
    ggplot2::geom_tile(colour = "gray92", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation,
                                     size   = abs(correlation))) +
    ggplot2::geom_point(
      data = dplyr::filter(data, significant),
      shape = 0, size = 10, stroke = 0.8, color = "black"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(guide = "none", range = c(1, 8),
                                   limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E", midpoint = 0,
      limits = c(-ran, ran), name = "Correlation", na.value = "gray95",
      guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(Biodomain_label), rows = ggplot2::vars(group),
      scales = "free", space = "free", switch = "y"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.x       = ggplot2::element_text(size = 12, face = "bold"),
      strip.text.y.left  = ggplot2::element_text(angle = 0, size = 12, face = "bold"),
      strip.background   = ggplot2::element_rect(fill = "gray95", color = "gray80"),
      axis.text.x        = ggplot2::element_text(angle = 90, hjust = 0, size = 11),
      axis.text.y        = ggplot2::element_text(size = 10),
      panel.spacing      = ggplot2::unit(0.1, "lines"),
      panel.border       = ggplot2::element_rect(color = "gray80", fill = NA),
      panel.grid         = ggplot2::element_blank(),
      legend.position    = "right"
    )
}


# ==============================================================================
# 14. ORA WITH BIODOMAIN ANNOTATION
# ==============================================================================

#' BioDomain-Aware GO Over-Representation Analysis
#'
#' @param dat      Data frame with symbol, log2FoldChange, group (padj-filtered)
#' @param set1     Character vector of group names to analyze
#' @param bio_ref  BioDomain mouse reference (annotated_biodomains_mmus.rds)
#' @param lab_ref  Domain color reference (domain_colors.csv)
#' @param universe Character vector of background gene symbols
#' @param org_db   OrgDb object (default: org.Mm.eg.db)
#' @return Data frame of enriched GO terms annotated with BioDomains
enrichORA_BD <- function(dat, set1,
                          bio_ref  = biodom_mouse,
                          lab_ref  = dom.lab,
                          universe = univ,
                          org_db   = org.Mm.eg.db) {

  process_group <- function(current_group) {
    message("ORA for group: ", current_group)

    group_dat <- dat %>%
      dplyr::filter(group == current_group, !is.na(symbol)) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)

    run_go <- function(genes, dir) {
      if (length(genes) == 0) return(data.frame())
      enr <- tryCatch(
        clusterProfiler::enrichGO(
          gene = genes, universe = universe, OrgDb = org_db,
          keyType = "SYMBOL", ont = "all", pAdjustMethod = "BH"
        ),
        error = function(e) NULL
      )
      res <- if (!is.null(enr)) as.data.frame(enr) else data.frame()
      if (nrow(res) > 0) res$dir <- dir
      return(res)
    }

    enr_ora <- dplyr::bind_rows(
      run_go(dplyr::filter(group_dat, log2FoldChange > 0)$symbol, "up"),
      run_go(dplyr::filter(group_dat, log2FoldChange < 0)$symbol, "dn")
    )
    if (nrow(enr_ora) == 0) return(NULL)

    bd_counts <- enr_ora %>%
      dplyr::left_join(bio_ref %>%
                         dplyr::select(Biodomain, Subdomain, ID = GO_ID),
                       by = "ID") %>%
      dplyr::mutate(
        Biodomain = tidyr::replace_na(Biodomain, "none"),
        Subdomain = dplyr::coalesce(Subdomain, Biodomain),
        model     = current_group
      ) %>%
      dplyr::group_by(Biodomain) %>%
      dplyr::summarise(n_sig_term = dplyr::n_distinct(ID), .groups = "drop")

    enr_ora %>%
      dplyr::left_join(bio_ref %>%
                         dplyr::select(Biodomain, Subdomain, ID = GO_ID),
                       by = "ID") %>%
      dplyr::mutate(
        Biodomain = tidyr::replace_na(Biodomain, "none"),
        Subdomain = dplyr::coalesce(Subdomain, Biodomain),
        model     = current_group
      ) %>%
      dplyr::left_join(lab_ref, by = c("Biodomain" = "domain")) %>%
      dplyr::left_join(bd_counts, by = "Biodomain") %>%
      dplyr::mutate(
        signed_logP = -log10(p.adjust),
        signed_logP = dplyr::if_else(dir == "dn", -1 * signed_logP, signed_logP),
        Biodomain   = forcats::fct_reorder(Biodomain, n_sig_term)
      ) %>%
      dplyr::arrange(Biodomain, p.adjust)
  }

  set1 %>% purrr::map_dfr(~process_group(.x))
}


#' Run Enrichment Pipeline for All DEG Groups
#'
#' @param df       DESeq2 results data frame (requires padj, symbol, log2FoldChange, group)
#' @param universe Character vector of background gene symbols
#' @param pval_cutoff Significance threshold for filtering (default: 0.05)
run_enrichment_pipeline <- function(df, universe, pval_cutoff = 0.05) {
  dat <- df %>%
    dplyr::filter(padj < pval_cutoff) %>%
    dplyr::select(symbol, EntrezGene, log2FoldChange, group) %>%
    tidyr::drop_na()

  enrichORA_BD(dat, unique(dat$group), universe = universe)
}


#' BioDomain ORA Violin + Jitter Plot
#'
#' @param enr.ora Output from enrichORA_BD()
ORA_plot_fn <- function(enr.ora) {
  ggplot2::ggplot(enr.ora, ggplot2::aes(x = signed_logP, y = model)) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.6,
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_violin(data = subset(enr.ora, dir == "up"),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = "grey90") +
    ggplot2::geom_violin(data = subset(enr.ora, dir == "dn"),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = "grey90") +
    ggplot2::geom_jitter(ggplot2::aes(size = Count, fill = color),
                         color = "black", shape = 21, alpha = 0.5,
                         stroke = 0.2, width = 0, height = 0.2) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(2, 7)) +
    ggplot2::labs(title = "ORA Results Categorized by Biodomains",
                  subtitle = "Directional Enrichment Analysis",
                  x = "Signed -log10(p-value)", y = "", size = "Gene Count") +
    ggplot2::facet_wrap(~Biodomain, scales = "free_x") +
    ggplot2::theme(
      axis.text.y      = ggplot2::element_text(size = 11, color = "black"),
      axis.text.x      = ggplot2::element_text(size = 11, color = "black"),
      axis.title       = ggplot2::element_text(size = 13, face = "bold"),
      plot.title       = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black"),
      strip.text       = ggplot2::element_text(size = 11, face = "bold"),
      legend.position  = "right"
    )
}


#' BioDomain Sub-Domain ORA Distribution Plot
#'
#' @param data  ORA data with Biodomain, Subdomain, model, color columns
#' @param title Plot title
plot_biodomain_distribution <- function(data,
                                         title = "BioDomain ORA Distribution") {
  ggplot2::ggplot(data, ggplot2::aes(x = signed_logP, y = Subdomain)) +
    ggplot2::geom_vline(xintercept = 0, lwd = 0.5, lty = 2, color = "gray50") +
    ggplot2::geom_violin(data = dplyr::filter(data, dir == "up"),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = NA) +
    ggplot2::geom_violin(data = dplyr::filter(data, dir == "dn"),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = NA) +
    ggplot2::geom_jitter(ggplot2::aes(size = Count, fill = color),
                         color = "grey20", shape = 21, alpha = 0.4,
                         width = 0, height = 0.25) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    ggplot2::facet_grid(rows = ggplot2::vars(Biodomain),
                        cols = ggplot2::vars(model),
                        scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title,
                  x = expression(paste(signed, " -log"[10], "(p-adj)")),
                  y = NULL) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "gray90", color = "gray70"),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 11, face = "bold"),
      strip.text.x      = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y       = ggplot2::element_text(size = 10),
      axis.text.x       = ggplot2::element_text(size = 11),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "bottom"
    )
}


# ==============================================================================
# 15. GSEA WITH BIODOMAIN ANNOTATION
# ==============================================================================

#' BioDomain-Aware GSEA
#'
#' @param dat     DESeq2 data frame with symbol, log2FoldChange, group
#' @param set1    Character vector of group names
#' @param bio_ref BioDomain mouse reference
#' @param lab_ref Domain color reference
#' @param org_db  OrgDb object (default: org.Mm.eg.db)
#' @return Combined data frame of GSEA results annotated with BioDomains
enrichGSEA_BD <- function(dat, set1,
                           bio_ref = biodom_mouse,
                           lab_ref = dom.lab,
                           org_db  = org.Mm.eg.db) {

  process_gsea <- function(current_group) {
    message("GSEA for group: ", current_group)

    ranked_list <- dat %>%
      dplyr::filter(group == current_group,
                    !is.na(symbol), !is.na(log2FoldChange)) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
      dplyr::distinct(symbol, .keep_all = TRUE) %>%
      { stats::setNames(.$log2FoldChange, .$symbol) }

    enr <- tryCatch(
      clusterProfiler::gseGO(
        geneList     = ranked_list, ont = "all", OrgDb = org_db,
        keyType      = "SYMBOL", pvalueCutoff = 1, seed = TRUE, verbose = FALSE
      ),
      error = function(e) {
        message("GSEA failed for ", current_group, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)

    enr_res <- as.data.frame(enr) %>%
      dplyr::left_join(bio_ref %>%
                         dplyr::select(Biodomain, Subdomain, ID = GO_ID),
                       by = "ID") %>%
      dplyr::mutate(
        Biodomain = tidyr::replace_na(Biodomain, "none"),
        Subdomain = dplyr::coalesce(Subdomain, Biodomain),
        model     = current_group
      )

    bd_counts <- enr_res %>%
      dplyr::group_by(Biodomain) %>%
      dplyr::summarise(n_sig_term = dplyr::n_distinct(ID), .groups = "drop")

    enr_res %>%
      dplyr::left_join(lab_ref, by = c("Biodomain" = "domain")) %>%
      dplyr::left_join(bd_counts, by = "Biodomain") %>%
      dplyr::mutate(Biodomain = forcats::fct_reorder(Biodomain, n_sig_term)) %>%
      dplyr::arrange(Biodomain, p.adjust)
  }

  set1 %>% purrr::map_dfr(~process_gsea(.x))
}


#' GSEA BioDomain Distribution Violin + Jitter Plot
#'
#' @param data Output from enrichGSEA_BD()
#' @param ran  NES axis limits (default: 3)
plot_gsea_biodomain_distribution <- function(data, ran = 3) {

  ggplot2::ggplot(data, ggplot2::aes(x = NES, y = model)) +
    ggplot2::geom_vline(xintercept = 0, lwd = 0.5, lty = 2, color = "gray50") +
    ggplot2::geom_violin(data = dplyr::filter(data, NES > 0),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = NA) +
    ggplot2::geom_violin(data = dplyr::filter(data, NES < 0),
                         ggplot2::aes(color = color), scale = "width",
                         alpha = 0.1, fill = NA) +
    ggplot2::geom_jitter(ggplot2::aes(size = -log10(p.adjust), fill = color),
                         color = "grey20", shape = 21, alpha = 0.4,
                         width = 0, height = 0.25) +
    ggplot2::scale_x_continuous(limits = c(-ran, ran),
                                 breaks = seq(-ran, ran, 1)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(1, 6),
                                   name = expression(paste(
                                     "-log"[10], "(p-adj)"))) +
    ggplot2::facet_wrap(~label) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "gray95", color = "gray80"),
      strip.text        = ggplot2::element_text(face = "bold", size = 10),
      axis.text         = ggplot2::element_text(size = 9),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "bottom"
    )
}


# ==============================================================================
# 16. VISUALIZATION: DEG BARPLOT
# ==============================================================================

#' DEG Directional Barplot (Up vs Down, Faceted by Sex)
#'
#' Requires `my_colors` to be defined in the calling environment:
#'   my_colors <- c("Up_DEGs" = "#D55E00", "Down_DEGs" = "#0072B2")
#'
#' @param data Data frame with columns: comparison, DEGs, Group, Sex
barplot_fn <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(x = comparison, y = DEGs, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge",
                      color = "black", linewidth = 0.3) +
    ggplot2::facet_wrap(~Sex, scales = "free_y") +
    ggplot2::scale_fill_manual(values = my_colors) +
    ggplot2::labs(
      title    = "Differentially Expressed Genes (DEGs)",
      subtitle = "FDR < 0.05 & |Log2FC| > 0",
      x = "", y = "Number of DEGs", fill = "Direction"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                                face = "bold", size = 11,
                                                color = "black"),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 11,
                                                color = "black"),
      axis.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                size = 16),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 12,
                                                face = "italic"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text       = ggplot2::element_text(face = "bold", size = 13),
      legend.position  = "top",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

# BioDomain reference — loaded in notebook, stored here for reference
human_biodomains <- tibble::tibble(
  Biodomain = c(
    "Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair",
    "Endolysosome", "Epigenetic", "Immune Response", "Lipid Metabolism",
    "Metal Binding and Homeostasis", "Mitochondrial Metabolism",
    "Myelination", "Oxidative Stress", "Proteostasis", "RNA Spliceosome",
    "Structural Stabilization", "Synapse", "Tau Homeostasis", "Vasculature"
  ),
  Biodomain_label = c(
    "Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair",
    "Endolysosome", "Epigenetic", "Immune\n Response", "Lipid Metabolism",
    "Metal Binding\n and Homeostasis", "Mitochondrial\n Metabolism",
    "Myelination", "Oxidative Stress", "Proteostasis", "RNA\n Spliceosome",
    "Structural Stabilization", "Synapse", "Tau\n Homeostasis", "Vasculature"
  )
)

# End of Helper_Functions.R
