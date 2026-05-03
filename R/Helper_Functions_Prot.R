# ==============================================================================
# Helper_Functions_Prot.R
# DIA Proteomics Analysis — Core Utility Functions
#
# Description:
#   Shared functions for proteomics QC, differential abundance analysis,
#   enrichment analysis, and visualization. Designed to work with any
#   case-control DIA proteomics study with sex and age stratification.
#
# Usage:
#   source("R/Helper_Functions_Prot.R")
#   All user-configurable settings (paths, organism, thresholds) should be
#   defined in the USER CONFIGURATION block in the main RMD notebook.
#
# Changelog:
#   - Added auto-detection of protein column start in anova_DEG_SampleType()
#   - Added safe_upset() with full empty-data guards
#   - Added map_symbols_to_entrez() utility for reliable gene ID mapping
#   - Updated kegg.fn() with NULL guards and NA Entrez filtering
#   - Added build_comparisons() to filter invalid group pairs automatically
#   - Added run_anova_all_comparisons() to move ANOVA loop out of notebook
#
# Author: Ravi S Pandey
# ==============================================================================


# ==============================================================================
# 1. SETUP & DEPENDENCIES
# ==============================================================================

needed.packages <- c(
  "ComplexHeatmap", "AnnotationDbi", "org.Mm.eg.db","circlize", "corrplot", "cowplot",
  "clusterProfiler", "dplyr", "data.table", "EnhancedVolcano",
  "enrichplot", "forcats", "ggplot2", "gridExtra", "ggpubr", "ggnewscale",
  "ggplotify", "ggrepel", "Hmisc", "janitor", "purrr",
  "RColorBrewer", "stringr", "scales", "tidyr", "tibble", "UpSetR",
  "gt", "gprofiler2", "readxl", "readr", "ggh4x", "pcaMethods"
)

# Load libraries silently — missing packages produce a clear error message
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
# 2. REFERENCE DATA LOADERS
# ==============================================================================

#' Load BioDomain Reference Data
#'
#' @param biodom_path        Path to annotated_biodomains .rds file
#' @param domain_colors_path Path to domain_colors .csv file
#' @return List with $biodom and $dom_lab
load_biodomain_refs <- function(biodom_path, domain_colors_path) {
  if (!file.exists(biodom_path))
    stop(paste("BioDomain RDS not found at:", biodom_path))
  if (!file.exists(domain_colors_path))
    stop(paste("Domain colors CSV not found at:", domain_colors_path))
  list(
    biodom  = readRDS(biodom_path),
    dom_lab = readr::read_csv(domain_colors_path, show_col_types = FALSE)
  )
}

#' Load Co-expression Module Reference Data
#'
#' @param module_rdata_path Path to .RData file containing module logFC objects
#' @return Invisibly; loads objects into the calling environment
load_module_refs <- function(module_rdata_path) {
  if (!file.exists(module_rdata_path))
    stop(paste("Module RData not found at:", module_rdata_path))
  load(module_rdata_path, envir = parent.frame())
}

# Canonical AD biological domain order (used for factor ordering in plots)
domain_order <- c(
  "Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair",
  "Endolysosome", "Epigenetic", "Immune Response", "Lipid Metabolism",
  "Metal Binding and Homeostasis", "Mitochondrial Metabolism",
  "Myelination", "Oxidative Stress", "Proteostasis", "RNA Spliceosome",
  "Structural Stabilization", "Synapse", "Tau Homeostasis", "Vasculature"
)

human_biodomains <- tibble::tibble(
  Biodomain = domain_order,
  Biodomain_label = c(
    "Apoptosis", "APP Metabolism", "Autophagy", "Cell Cycle", "DNA Repair",
    "Endolysosome", "Epigenetic", "Immune\n Response", "Lipid Metabolism",
    "Metal Binding\n and Homeostasis", "Mitochondrial\n Metabolism",
    "Myelination", "Oxidative Stress", "Proteostasis", "RNA\n Spliceosome",
    "Structural Stabilization", "Synapse", "Tau\n Homeostasis", "Vasculature"
  )
)


# ==============================================================================
# 3. CORE DIFFERENTIAL ANALYSIS FUNCTIONS
# ==============================================================================

#' Build Valid Pairwise Comparisons from Available Groups
#'
#' Generates all case-vs-control comparisons stratified by sex and age,
#' then filters to only pairs where BOTH groups exist in the data.
#' This prevents ANOVA errors from being run on empty subsets.
#'
#' @param df               Data frame containing a `Group` column
#' @param control_genotype String. Name of the control genotype
#' @param case_genotypes   Character vector of case genotype names
#' @param sex_levels       Character vector of sex levels (e.g., c("Male","Female"))
#' @param age_groups       Numeric vector of age timepoints (e.g., c(4, 12))
#' @return Data frame with columns: control, case (only valid pairs)
build_comparisons <- function(df, control_genotype, case_genotypes,
                              sex_levels, age_groups) {
  available_groups <- unique(df$Group)
  
  comparisons <- expand.grid(
    case_genotype = case_genotypes,
    sex           = sex_levels,
    age           = age_groups,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      control = paste0(control_genotype, "_", sex, "_", age, "M"),
      case    = paste0(case_genotype,    "_", sex, "_", age, "M")
    ) %>%
    # Only keep pairs where both groups actually exist in the data
    dplyr::filter(
      control %in% available_groups,
      case    %in% available_groups
    ) %>%
    dplyr::select(control, case)
  
  if (nrow(comparisons) == 0)
    warning("No valid comparisons found. Check that group names match your data exactly.")
  
  message("Valid comparisons found: ", nrow(comparisons))
  return(comparisons)
}


#' ANOVA + Tukey HSD for Differential Abundance
#'
#' For each protein, runs a one-way ANOVA across SampleType groups and
#' extracts Tukey HSD post-hoc p-values and fold-changes.
#' Protein columns are auto-detected as the first numeric columns,
#' so no hardcoded column index is needed.
#'
#' @param data             Data frame with a `SampleType` column and protein columns
#' @param protein_col_start Integer or NULL. Column index where proteins begin.
#'                         If NULL (default), auto-detected from first numeric column.
#' @return Data frame with F-value, raw p-value, Tukey adjusted p-value, and diff per protein
anova_DEG_SampleType <- function(data, protein_col_start = NULL) {
  
  # Auto-detect protein column start from first numeric column if not specified
  if (is.null(protein_col_start)) {
    protein_col_start <- min(which(sapply(data, is.numeric)))
    message("Auto-detected protein columns starting at column: ", protein_col_start)
  }
  
  # Remove proteins where any SampleType group is entirely NA
  cols_to_remove <- data %>%
    dplyr::group_by(SampleType) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~all(is.na(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-SampleType) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~any(.x))) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "col", values_to = "is_na") %>%
    dplyr::filter(is_na) %>%
    dplyr::pull(col)
  
  if (length(cols_to_remove) > 0)
    message("Removing ", length(cols_to_remove), " proteins with complete NA in at least one group.")
  
  df_cleaned <- data %>% dplyr::select(-dplyr::any_of(cols_to_remove))
  
  # Run ANOVA + Tukey for each protein column
  res_list <- lapply(protein_col_start:ncol(df_cleaned), function(i) {
    tryCatch({
      aov_obj <- aov(df_cleaned[[i]] ~ SampleType, data = df_cleaned)
      a_res   <- anova(aov_obj)
      t_res   <- as.data.frame(TukeyHSD(aov_obj)$SampleType)
      
      c(Protein = colnames(df_cleaned)[i],
        F_Val   = a_res$F[1],
        P_Val   = a_res$Pr[1],
        p_adj   = t_res$`p adj`,
        diff    = t_res$diff)
    }, error = function(e) NULL)
  })
  
  # Remove NULL results (proteins that failed)
  res_list  <- Filter(Negate(is.null), res_list)
  final_df  <- as.data.frame(do.call(rbind, res_list))
  rownames(final_df) <- final_df$Protein
  return(final_df[, -1])
}


#' Run ANOVA Across All Valid Comparisons
#'
#' Wraps the ANOVA loop that would otherwise live in the notebook.
#' Takes a pre-built comparisons table (from build_comparisons()) and
#' returns both a list and a combined data frame of results.
#'
#' @param df          Data frame with Group column and protein abundance columns
#' @param comparisons Data frame with columns: control, case (from build_comparisons())
#' @return List with $results_list (named list) and $results_df (combined data frame)
run_anova_all_comparisons <- function(df, comparisons) {
  
  if (nrow(comparisons) == 0)
    stop("No valid comparisons to run. Check build_comparisons() output.")
  
  DE_results.list <- list()
  DE_results.df   <- data.frame()
  
  for (i in 1:nrow(comparisons)) {
    ctrl       <- comparisons[i, "control"]
    case       <- comparisons[i, "case"]
    model_name <- paste0(case, "_vs_", ctrl)
    
    message("Running ANOVA: ", model_name)
    
    data_subset <- df %>%
      dplyr::filter(Group %in% c(ctrl, case)) %>%
      dplyr::rename(SampleType = Group)
    
    # Skip if either group has fewer than 3 samples
    group_counts <- table(droplevels(data_subset$SampleType))
    if (any(group_counts < 2)) {
      message("  Skipping ", model_name, " — insufficient samples (n < 2) in at least one group.")
      next
    }
    
    anova_results <- tryCatch({
      anova_DEG_SampleType(data_subset) %>%
        as.data.frame() %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
        tibble::rownames_to_column(var = "Symbol")
    }, error = function(e) {
      message("  ANOVA failed for ", model_name, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(anova_results)) next
    
    colnames(anova_results) <- c("Symbol", "F.Value", "Pr.F", "padj", "diff")
    
    DE_results.df                   <- rbind(DE_results.df,
                                             anova_results %>% dplyr::mutate(model = model_name))
    DE_results.list[[model_name]]   <- anova_results
  }
  
  message("Completed: ", length(DE_results.list), " of ", nrow(comparisons), " comparisons.")
  return(list(results_list = DE_results.list, results_df = DE_results.df))
}


# ==============================================================================
# 4. GENE ID MAPPING UTILITIES
# ==============================================================================

#' Map Gene Symbols to Entrez IDs
#'
#' Resolves gene symbols to Entrez IDs using AnnotationDbi.
#' Operates on unique symbols only (faster), then maps back to the full vector.
#' Always call this OUTSIDE of mutate() to avoid scoping issues.
#'
#' @param symbols  Character vector of gene symbols
#' @param org_db   OrgDb object, already loaded (e.g., org.Mm.eg.db)
#' @return Named character vector of Entrez IDs (NA where mapping fails)
#'
#' @examples
#' # In your notebook USER CONFIGURATION block:
#' library(organism_db, character.only = TRUE)
#' org_db_obj <- get(organism_db)
#'
#' # Then use:
#' entrez_map <- map_symbols_to_entrez(female_degs$Symbol, org_db_obj)
#' dat_f <- female_degs %>% mutate(EntrezGene = entrez_map[Symbol])
map_symbols_to_entrez <- function(symbols, org_db) {
  if (is.character(org_db))
    stop("Pass the OrgDb object directly, not a string. ",
         "Use: org_db_obj <- get(organism_db) first.")
  
  AnnotationDbi::mapIds(
    org_db,
    keys      = unique(as.character(symbols)),
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
}


#' Add Entrez IDs to a DEG Data Frame
#'
#' Convenience wrapper around map_symbols_to_entrez() that adds
#' an EntrezGene column directly to a data frame.
#'
#' @param deg_df  Data frame with a Symbol column
#' @param org_db  OrgDb object (already loaded)
#' @return deg_df with an added EntrezGene column
add_entrez_ids <- function(deg_df, org_db) {
  entrez_map <- map_symbols_to_entrez(deg_df$Symbol, org_db)
  deg_df %>% dplyr::mutate(EntrezGene = entrez_map[Symbol])
}


# ==============================================================================
# 5. OVERLAP & UPSET ANALYSIS
# ==============================================================================

#' Extract Up/Down Regulated Protein Lists Per Model
#'
#' Generates named lists of up- and down-regulated proteins for each model,
#' suitable for UpSetR visualization. Assigns all.up and all.down
#' to the global environment for compatibility with UpSetR.
#'
#' @param gp1      Character vector of model/comparison names
#' @param deg_data Data frame with columns: Symbol, model, diff
overlap.fn <- function(gp1, deg_data) {
  dat.up <- purrr::map(gp1, ~ deg_data %>%
                         dplyr::filter(model == .x, diff > 0) %>%
                         dplyr::pull(Symbol))
  dat.dn <- purrr::map(gp1, ~ deg_data %>%
                         dplyr::filter(model == .x, diff < 0) %>%
                         dplyr::pull(Symbol))
  names(dat.up) <- names(dat.dn) <- gp1
  all.up  <<- dat.up
  all.down <<- dat.dn
}


#' Safe UpSet Plot with Empty Data Guards
#'
#' Wraps UpSetR::upset() with checks for empty gene lists, insufficient groups,
#' and missing significance hits — all common with small or sparse datasets.
#'
#' @param deg_data   Significance-filtered DEG data frame (padj already applied upstream)
#' @param direction  "up" or "down" — which regulated direction to plot
#' @param bar_color  Bar color for the upset plot
#' @param x_label    X axis label
#' @param min_groups Minimum number of non-empty groups required to draw the plot (default: 2)
#' @return Invisibly NULL if skipped, otherwise renders the upset plot
safe_upset <- function(deg_data,
                       direction  = "up",
                       bar_color  = "#164B6E",
                       x_label    = "Total Proteins",
                       min_groups = 2) {
  
  # Guard 1: no significant proteins at all
  if (is.null(deg_data) || nrow(deg_data) == 0) {
    message("Upset plot skipped: no significant proteins found.")
    return(invisible(NULL))
  }
  
  models <- unique(deg_data$model)
  
  # Guard 2: only one model — nothing to overlap
  if (length(models) < min_groups) {
    message("Upset plot skipped: need at least ", min_groups, " models, found ", length(models), ".")
    return(invisible(NULL))
  }
  
  overlap.fn(models, deg_data)
  plot_list <- if (direction == "up") all.up else all.down
  
  # Guard 3: all lists are empty
  non_empty <- sapply(plot_list, length) > 0
  if (sum(non_empty) == 0) {
    message("Upset plot skipped: no ", direction, "regulated proteins in any group.")
    return(invisible(NULL))
  }
  
  # Guard 4: only one non-empty group — no overlap possible
  if (sum(non_empty) < min_groups) {
    message("Upset plot skipped: need at least ", min_groups,
            " groups with proteins, found ", sum(non_empty), ".")
    return(invisible(NULL))
  }
  
  UpSetR::upset(
    UpSetR::fromList(plot_list),
    nsets          = length(models),
    main.bar.color = bar_color,
    sets.x.label   = x_label,
    order.by       = "freq"
  )
}


# ==============================================================================
# 6. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Differential Abundance Barplot
#'
#' @param data        Data frame with columns: comparison, DEPs, Direction, Sex
#' @param pval_cutoff Numeric. Used in subtitle label
#' @param lfc_cutoff  Numeric. Used in subtitle label
barplot_fn <- function(data, pval_cutoff = 0.05, lfc_cutoff = 0) {
  my_colors <- c("Up" = "#D55E00", "Down" = "#0072B2")
  
  ggplot2::ggplot(data, ggplot2::aes(x = comparison, y = DEPs, fill = Direction)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge",
                      color = "black", linewidth = 0.3) +
    ggplot2::facet_wrap(~Sex, scales = "free_y") +
    ggplot2::scale_fill_manual(values = my_colors) +
    ggplot2::labs(
      title    = "Differentially Abundant Proteins (DAPs)",
      subtitle = paste0("adjP < ", pval_cutoff, " | |Log2FC| > ", lfc_cutoff),
      x = "", y = "Number of DAPs", fill = "Direction"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                               face = "bold", size = 11, color = "black"),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 11, color = "black"),
      axis.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 12, face = "italic"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text       = ggplot2::element_text(face = "bold", size = 13),
      legend.position  = "top",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}


#' Protein Fold-Change Expression Bar Plot
#'
#' Bar chart of log2 fold-change for selected proteins across comparisons.
#' Asterisks mark entries below the p-value cutoff.
#'
#' @param deg_df      DEG results filtered to one sex
#' @param proteins    Character vector of protein symbols to plot
#' @param sex_label   String used in plot title (e.g., "Females")
#' @param pval_cutoff Numeric. Significance threshold for asterisk display
plot_protein_expression <- function(deg_df, proteins, sex_label, pval_cutoff = 0.05) {
  
  plot_data <- deg_df %>%
    dplyr::filter(Symbol %in% proteins) %>%
    dplyr::mutate(significant = padj < pval_cutoff) %>%
    dplyr::filter(!is.na(model))
  
  if (nrow(plot_data) == 0) {
    message("No data found for proteins: ", paste(proteins, collapse = ", "))
    return(invisible(NULL))
  }
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = model, y = diff, fill = model)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(significant, "*", "")),
      vjust = ifelse(plot_data$diff > 0, -0.2, 1.2),
      size = 8, fontface = "bold"
    ) +
    ggplot2::facet_wrap(~Symbol, scales = "free_y") +
    ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    ggplot2::labs(
      title    = paste("Protein Abundance Changes:", sex_label),
      subtitle = paste0("* = adjP < ", pval_cutoff),
      x = "Comparison", y = "log2 Fold Change", fill = "Comparison"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text       = ggplot2::element_text(face = "bold", size = 12),
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}


#' Raw Abundance Boxplot for Selected Proteins
#'
#' @param traits_df Data frame with Genotype, Sex, Age, and protein abundance columns
#' @param proteins  Character vector of protein symbols to plot
#' @param group_col Name of the grouping column (default: "Genotype")
plot_protein_raw_abundance <- function(traits_df, proteins, group_col = "Genotype") {
  
  df_long <- traits_df %>%
    dplyr::select(dplyr::any_of(c("animalName", group_col, "Sex", "Age")),
                  dplyr::all_of(proteins)) %>%
    tidyr::pivot_longer(cols     = dplyr::all_of(proteins),
                        names_to  = "Protein",
                        values_to = "LogAbundance")
  
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = .data[[group_col]], y = LogAbundance,
                               fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.colour = "black", outlier.shape = 16,
                          outlier.size = 2, notch = FALSE) +
    ggplot2::ylab("log2(abundance)") +
    ggplot2::facet_wrap(~Sex + Protein, scales = "free") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(face = "bold", size = 12,
                                               angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "gray95"),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "none"
    ) +
    ggplot2::scale_fill_brewer(palette = "Spectral")
}


# ==============================================================================
# 7. CO-EXPRESSION MODULE CORRELATION ANALYSIS
# ==============================================================================

#' Correlate Mouse DEG Results with Human AD Co-expression Modules
#'
#' @param deg_df    DEG results (Symbol, model, diff)
#' @param human_ref Human reference with columns: Symbol, Module, AD.Control
#' @return Data frame of Pearson correlations per Module x Model
run_module_correlation <- function(deg_df, human_ref) {
  deg_df %>%
    dplyr::select(Symbol, model, diff) %>%
    dplyr::inner_join(human_ref, by = "Symbol") %>%
    dplyr::group_by(Module, model) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    tidyr::nest(data = c(Symbol, diff, AD.Control)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["AD.Control"]],
                                                method = "pearson")),
      correlation = purrr::map_dbl(cor_test, "estimate"),
      p_value     = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      FDR         = p.adjust(p_value, method = "fdr"),
      significant = p_value < 0.05
    )
}


#' Module Correlation Heat-Dot Plot
#'
#' @param cor_data    Output from run_module_correlation()
#' @param ran         Color scale limit (auto-computed if NULL)
#' @param order_model Character vector defining row order of models
#' @param module_order Character vector defining column order of modules
#' @param title       String. Plot title
plot_module_corrplot <- function(cor_data, ran = NULL, order_model = NULL,
                                 module_order = NULL, title = "Human AD Modules") {
  if (is.null(ran))
    ran <- ceiling(max(abs(range(cor_data$correlation, na.rm = TRUE))) * 100) / 100
  
  plot_df <- cor_data
  if (!is.null(module_order))
    plot_df <- plot_df %>% dplyr::mutate(Module = factor(Module, levels = module_order))
  if (!is.null(order_model))
    plot_df <- plot_df %>%
    dplyr::mutate(model = factor(model, levels = order_model),
                  model = forcats::fct_rev(model))
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Module, y = model)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation, size = abs(correlation))) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, significant),
                        color = "black", shape = 0, size = 10, stroke = 1.2) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E",
      midpoint = 0, limits = c(-ran, ran), name = "Correlation (r)",
      guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15)
    ) +
    ggplot2::labs(x = "Human Brain Proteomics Module", y = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 90, hjust = 0, size = 10, face = "bold"),
      axis.text.y     = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      panel.grid      = ggplot2::element_blank()
    )
}


#' Scatter Plot of Module Correlations Per Model
#'
#' @param deg_df       DEG results (Symbol, model, diff)
#' @param human_ref    Human reference module data
#' @param mods_to_plot Character vector of module names to visualize
#' @param sex_label    Label for y-axis
plot_module_scatter_correlations <- function(deg_df, human_ref,
                                             mods_to_plot, sex_label) {
  for (mod in mods_to_plot) {
    plot_data <- deg_df %>%
      dplyr::select(Symbol, model, diff) %>%
      dplyr::inner_join(human_ref, by = "Symbol") %>%
      dplyr::filter(Module %in% mod)
    
    if (nrow(plot_data) > 0) {
      p <- ggpubr::ggscatter(
        plot_data, x = "AD.Control", y = "diff",
        border = TRUE, label = "Symbol",
        add = "reg.line", conf.int = TRUE,
        xlab = paste("Human AD Module:", mod),
        ylab = paste(sex_label, "Mouse models (logFC)"),
        font.label = 12, repel = FALSE
      ) +
        ggplot2::theme_bw() +
        ggpubr::stat_cor(method = "pearson") +
        ggpubr::grids(linetype = "dashed") +
        ggplot2::facet_wrap(~model, scales = "free")
      print(p)
    } else {
      message("No data for module: ", mod)
    }
  }
}


# ==============================================================================
# 7.1 LINEAR MODELING UTILITIES
# ==============================================================================

#' Extract Estimates and P-values from LM Summary Objects
#'
#' @param lm_summaries Named list of summary(lm(...)) objects, one per protein
#' @param coef_names   Character vector of coefficient names to extract
#' @return Long-format data frame: Symbol, Variant, Estimate, P.Value
extract_lm_stats <- function(lm_summaries, coef_names) {
  estimates <- as.data.frame(
    dplyr::bind_rows(lapply(lm_summaries, function(x) x$coefficients[, "Estimate"]))
  )
  pvals <- as.data.frame(
    dplyr::bind_rows(lapply(lm_summaries, function(x) x$coefficients[, "Pr(>|t|)"]))
  )
  colnames(estimates) <- colnames(pvals) <- coef_names
  
  estimates %>%
    dplyr::mutate(Symbol = names(lm_summaries)) %>%
    dplyr::select(-`(Intercept)`) %>%
    tidyr::pivot_longer(cols = -Symbol, names_to = "Variant", values_to = "Estimate") %>%
    dplyr::left_join(
      pvals %>%
        dplyr::mutate(Symbol = names(lm_summaries)) %>%
        dplyr::select(-`(Intercept)`) %>%
        tidyr::pivot_longer(cols = -Symbol, names_to = "Variant", values_to = "P.Value"),
      by = c("Symbol", "Variant")
    )
}


#' Correlate LM Estimates with Human Module Signatures
#'
#' @param lm_results      Output from extract_lm_stats()
#' @param human_ref       Human reference with Symbol, Module, AD.Control
#' @param ordered_variant Character vector defining factor order of Variant
#' @param module_order    Character vector defining factor order of Module
#' @return Data frame with correlation and significance per Module x Variant
correlate_lm_with_human_modules <- function(lm_results, human_ref,
                                            ordered_variant = NULL,
                                            module_order = NULL) {
  cor_results <- lm_results %>%
    dplyr::select(-P.Value) %>%
    dplyr::inner_join(human_ref, by = "Symbol") %>%
    dplyr::group_by(Module, Variant) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    tidyr::nest(data = c(Symbol, Estimate, AD.Control)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["AD.Control"]],
                                                method = "pearson")),
      estimate    = purrr::map_dbl(cor_test, "estimate"),
      p_value     = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      padj         = p.adjust(p_value, method = "fdr"),
      significant  = padj < 0.05,
      significant1 = p_value < 0.05
    ) %>%
    dplyr::select(Module, Variant, correlation = estimate,
                  p_value, padj, significant, significant1)
  
  if (!is.null(module_order))
    cor_results <- cor_results %>%
      dplyr::mutate(Module = factor(Module, levels = module_order))
  if (!is.null(ordered_variant))
    cor_results <- cor_results %>%
      dplyr::mutate(Variant = factor(Variant, levels = ordered_variant),
                    Variant = forcats::fct_rev(Variant))
  cor_results
}


# ==============================================================================
# 8. BIODOMAIN ANALYSIS
# ==============================================================================

#' Correlate Model DEG Fold-Changes with Human BioDomain Reference
#'
#' @param deg_fc      DEG results (Symbol, model, diff)
#' @param human_ref   Human reference (Gene, Term, ampad_fc)
#' @param sig         p-value cutoff for significance flag
#' @param order_model Character vector for model factor ordering
corr_function_biodomain_protein <- function(deg_fc, human_ref, sig = 0.05,
                                            order_model = NULL) {
  dat <- deg_fc %>%
    dplyr::rename(Gene = Symbol) %>%
    dplyr::select(Gene, model, diff) %>%
    dplyr::inner_join(human_ref, by = "Gene") %>%
    dplyr::select(Term, model, Gene, diff, ampad_fc) %>%
    dplyr::group_by(Term, model) %>%
    tidyr::nest(data = c(Gene, diff, ampad_fc)) %>%
    dplyr::mutate(
      cor_test  = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["ampad_fc"]],
                                              method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(biodom = Term, model, correlation = estimate, p_value, significant) %>%
    dplyr::mutate(biodom = factor(biodom, levels = domain_order))
  
  if (!is.null(order_model))
    dat <- dat %>%
      dplyr::mutate(model = factor(model, levels = order_model),
                    model = forcats::fct_rev(model))
  dat
}


#' BioDomain Model Correlation Dot Plot
#'
#' @param data Output from corr_function_biodomain_protein()
biodom_model_corrplot <- function(data) {
  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100
  
  ggplot2::ggplot(data, ggplot2::aes(x = biodom, y = model)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation, size = abs(correlation)),
                        alpha = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(data, significant),
                        color = "black", shape = 0, size = 10, stroke = 1.1) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none", limits = c(0, ran)) +
    ggplot2::scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E",
      midpoint = 0, limits = c(-ran, ran), name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE, barheight = 10,
                                      frame.colour = "black")
    ) +
    ggplot2::facet_grid(cols = vars(ref), scales = "free", space = "free", switch = "y") +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "grey95", color = "black"),
      strip.text.x      = ggplot2::element_text(face = "bold", size = 12),
      strip.text.y.left = ggplot2::element_text(face = "bold", size = 12, angle = 0),
      axis.text.x       = ggplot2::element_text(angle = 90, hjust = 0, size = 11, color = "black"),
      axis.text.y       = ggplot2::element_text(size = 11, color = "black"),
      axis.ticks        = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "right"
    )
}


#' Correlate Sub-BioDomain with DEG Fold-Changes
#'
#' @param deg_fc      DEG results (Symbol, model, diff)
#' @param human_ref   Sub-biodomain reference (Gene, Term, ampad_fc, Biodomain)
#' @param sig         p-value cutoff
#' @param order_model Factor order for models
corr_function_subdomain_protein <- function(deg_fc, human_ref, sig = 0.05,
                                            order_model = NULL) {
  dat <- deg_fc %>%
    dplyr::rename(Gene = Symbol) %>%
    dplyr::select(Gene, model, diff) %>%
    dplyr::inner_join(human_ref, by = "Gene") %>%
    dplyr::select(Term, model, Gene, diff, ampad_fc, Biodomain) %>%
    dplyr::group_by(Term, model) %>%
    tidyr::drop_na() %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    tidyr::nest(data = c(Gene, diff, ampad_fc)) %>%
    dplyr::mutate(
      cor_test  = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["ampad_fc"]],
                                              method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(Biodomain, biodom = Term, model,
                  correlation = estimate, p_value, significant) %>%
    dplyr::mutate(Biodomain = factor(Biodomain, levels = domain_order))
  
  if (!is.null(order_model))
    dat <- dat %>%
      dplyr::mutate(model = factor(model, levels = order_model),
                    model = forcats::fct_rev(model))
  dat
}


#' Sub-BioDomain Portrait Correlation Plot
#'
#' @param data Output from corr_function_subdomain_protein()
plot_biodom_portrait_correlations <- function(data) {
  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100
  
  ggplot2::ggplot(data, ggplot2::aes(y = biodom, x = model)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation, size = abs(correlation))) +
    ggplot2::geom_point(data = dplyr::filter(data, significant),
                        color = "black", shape = 0, size = 9, stroke = 1.1) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran), range = c(1, 8)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation",
      guide = ggplot2::guide_colorbar(ticks = FALSE, title.position = "top",
                                      title.hjust = 0.5)
    ) +
    ggplot2::facet_grid(rows = vars(Biodomain), scales = "free",
                        space = "free", switch = "y") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 11, face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "gray96", color = "white"),
      axis.ticks        = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(angle = 90, hjust = 0,
                                                size = 11, color = "black"),
      axis.text.y       = ggplot2::element_text(size = 10, color = "black"),
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "right"
    )
}


#' Correlate LM Variants with BioDomain Reference
#'
#' @param data      LM results (Symbol, Variant, Estimate)
#' @param human_ref BioDomain reference (Gene, Term, ampad_fc)
#' @param sig       p-value cutoff
corr_function_lm_biodomain <- function(data, human_ref, sig = 0.05) {
  data %>%
    dplyr::rename(Gene = Symbol) %>%
    dplyr::select(Gene, Variant, Estimate) %>%
    dplyr::inner_join(human_ref, by = "Gene") %>%
    dplyr::group_by(Term, Variant) %>%
    tidyr::nest(data = c(Gene, Estimate, ampad_fc)) %>%
    dplyr::mutate(
      cor_test  = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["ampad_fc"]],
                                              method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(biodom = Term, Variant, correlation = estimate, p_value, significant) %>%
    dplyr::mutate(module = factor(biodom, levels = domain_order))
}


#' BioDomain Variant Correlation Dot Plot
#'
#' @param data Output from corr_function_lm_biodomain()
variant_corrplot_biodom <- function(data) {
  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100
  
  ggplot2::ggplot(data, ggplot2::aes(x = biodom, y = Variant)) +
    ggplot2::geom_tile(colour = "black", fill = "white") +
    ggplot2::geom_point(ggplot2::aes(colour = correlation, size = abs(correlation))) +
    ggplot2::geom_point(data = dplyr::filter(data, significant),
                        ggplot2::aes(colour = correlation),
                        color = "black", shape = 0, size = 9, stroke = 1.1) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(guide = "none", limits = c(0, ran), range = c(1, 8)) +
    ggplot2::scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation", guide = ggplot2::guide_colorbar(ticks = FALSE)
    ) +
    ggplot2::facet_grid(rows = vars(group), cols = vars(ref),
                        scales = "free", space = "free", switch = "y") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "grey95", color = "black"),
      strip.text.x      = ggplot2::element_text(angle = 0, size = 14, face = "bold"),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 14, face = "bold"),
      axis.ticks        = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(angle = 90, hjust = 0,
                                                size = 14, color = "black"),
      axis.text.y       = ggplot2::element_text(size = 14, color = "black"),
      panel.grid        = ggplot2::element_blank(),
      legend.position   = "right"
    )
}


#' Correlate LM Variants with Sub-BioDomain Reference
corr_function_lm_subbiodomain <- function(lm_results, human_ref, sig = 0.05) {
  lm_results %>%
    dplyr::rename(Gene = Symbol) %>%
    dplyr::select(Gene, Variant, Estimate) %>%
    dplyr::inner_join(human_ref, by = "Gene") %>%
    dplyr::group_by(Term, Variant) %>%
    tidyr::drop_na() %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    tidyr::nest(data = c(Gene, Estimate, ampad_fc)) %>%
    dplyr::mutate(
      cor_test  = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["ampad_fc"]],
                                              method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(human_biodomains, by = "Biodomain") %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(module = Biodomain, biodom = Term, Variant,
                  correlation = estimate, p_value, significant, Biodomain_label) %>%
    dplyr::mutate(module = factor(module, levels = domain_order))
}


# ==============================================================================
# 9. ENRICHMENT UTILITIES — ORA
# ==============================================================================

#' GO Over-Representation Analysis with BioDomain Mapping
#'
#' @param dat         DEG data frame (Symbol, diff, padj, model)
#' @param bio_ref     BioDomain reference (from load_biodomain_refs()$biodom)
#' @param lab_ref     Domain color labels (from load_biodomain_refs()$dom_lab)
#' @param org_db      OrgDb object (already loaded, NOT a string)
#' @param pval_cutoff Significance threshold for gene selection
#' @param universe    Optional background gene set; auto-generated if NULL
#' @return Data frame of enriched GO terms annotated with BioDomains
enrichORA_BD_Protein <- function(dat, bio_ref, lab_ref, org_db,
                                 pval_cutoff = 0.05, universe = NULL) {
  if (is.character(org_db))
    stop("Pass the OrgDb object directly. Use: org_db_obj <- get(organism_db) first.")
  
  if (is.null(universe)) {
    message("Generating background universe from OrgDb...")
    universe <- AnnotationDbi::mapIds(
      org_db,
      keys     = AnnotationDbi::keys(org_db, "SYMBOL"),
      column   = "SYMBOL",
      keytype  = "SYMBOL",
      multiVals = "first"
    ) %>% intersect(unique(dat$Symbol))
  }
  
  gene_lists <- dat %>%
    dplyr::select(symbol = Symbol, logFC = diff, padj, model) %>%
    dplyr::filter(!is.na(symbol), padj <= pval_cutoff) %>%
    dplyr::mutate(dir = ifelse(logFC > 0, "up", "dn")) %>%
    dplyr::group_by(model, dir) %>%
    dplyr::summarise(genelist = list(symbol), .groups = "drop")
  
  if (nrow(gene_lists) == 0) {
    message("No significant genes found for ORA.")
    return(NULL)
  }
  
  message("Running GO Enrichment (ORA)...")
  go_results <- gene_lists %>%
    dplyr::mutate(
      enr.result = purrr::map(genelist, function(genes) {
        res <- tryCatch({
          clusterProfiler::enrichGO(
            gene = genes, universe = universe, OrgDb = org_db,
            keyType = "SYMBOL", ont = "all", pAdjustMethod = "BH"
          )
        }, error = function(e) NULL)
        if (is.null(res)) return(data.frame())
        as.data.frame(res)
      })
    ) %>%
    tidyr::unnest(enr.result) %>%
    dplyr::select(-genelist)
  
  if (nrow(go_results) == 0) {
    message("No significant GO terms enriched.")
    return(NULL)
  }
  
  message("Mapping GO terms to BioDomains...")
  go_results %>%
    dplyr::inner_join(
      bio_ref %>% dplyr::select(ID = GO_ID, Biodomain, Subdomain),
      by = "ID", relationship = "many-to-many"
    ) %>%
    dplyr::left_join(lab_ref, by = c("Biodomain" = "domain")) %>%
    dplyr::mutate(
      Biodomain   = tidyr::replace_na(Biodomain, "none"),
      Subdomain   = tidyr::replace_na(Subdomain, "none"),
      signed_logP = -log10(p.adjust),
      signed_logP = ifelse(dir == "dn", -1 * signed_logP, signed_logP)
    ) %>%
    dplyr::group_by(Biodomain) %>%
    dplyr::mutate(n_sig = dplyr::n_distinct(ID)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(n_sig), p.adjust)
}


#' ORA BioDomain Violin + Jitter Plot
#'
#' @param enr.ora Output from enrichORA_BD_Protein()
enr.ORA.bd_plot <- function(enr.ora) {
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
    ggplot2::labs(
      title    = "ORA Results Categorized by Biodomains",
      subtitle = "Directional Enrichment Analysis",
      x = "Signed -log10(p-value)", y = "", size = "Gene Count"
    ) +
    ggplot2::facet_wrap(~Biodomain, scales = "free_x") +
    ggplot2::theme(
      axis.text.y      = ggplot2::element_text(size = 11, face = "bold", color = "black"),
      axis.text.x      = ggplot2::element_text(size = 11, color = "black"),
      axis.title       = ggplot2::element_text(size = 13, face = "bold"),
      plot.title       = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black"),
      strip.text       = ggplot2::element_text(size = 11, face = "bold"),
      legend.position  = "right"
    )
}


#' Sub-BioDomain ORA Distribution Plot
#'
#' @param data  ORA data with Biodomain, Subdomain, model columns
#' @param title Plot title
plot_subdomain_distribution <- function(data, title = "Sub-BioDomain ORA Distribution") {
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
    ggplot2::facet_grid(rows = vars(Biodomain), cols = vars(model),
                        scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title,
                  x = expression(paste(signed, " -log"[10], "(p-adj)")), y = NULL) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "gray90", color = "gray70"),
      strip.text.y.left = ggplot2::element_text(angle = 0, size = 11, face = "bold"),
      strip.text.x      = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y       = ggplot2::element_text(size = 10),
      axis.text.x       = ggplot2::element_text(size = 11),
      panel.border      = ggplot2::element_rect(color = "gray90", fill = NA),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "bottom"
    )
}


# ==============================================================================
# 10. GSEA WITH BIODOMAIN ANNOTATION
# ==============================================================================

#' Gene Set Enrichment Analysis with BioDomain Mapping
#'
#' @param dat     DEG data frame (Symbol, diff, model)
#' @param bio_ref BioDomain reference
#' @param lab_ref Domain color labels
#' @param org_db  OrgDb object (already loaded, NOT a string)
#' @param seed    Integer. Random seed for reproducibility
#' @return Data frame of GSEA results annotated with BioDomains
enrichGSEA_BD_Protein <- function(dat, bio_ref, lab_ref, org_db, seed = 1111) {
  if (is.character(org_db))
    stop("Pass the OrgDb object directly. Use: org_db_obj <- get(organism_db) first.")
  
  set.seed(seed)
  message("Running GSEA with GO gene sets...")
  
  enr_GSEA <- dat %>%
    dplyr::select(symbol = Symbol, log2FoldChange = diff, padj, model) %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      gl = list(setNames(log2FoldChange, symbol) %>% sort(decreasing = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      gse = purrr::map(gl, ~ clusterProfiler::gseGO(
        .x, ont = "all", OrgDb = org_db, keyType = "SYMBOL",
        pvalueCutoff = 1, eps = 0, seed = TRUE
      )),
      res = purrr::map(gse, ~as.data.frame(.x))
    )
  
  message("Mapping GSEA results to BioDomains...")
  enr_GSEA %>%
    dplyr::select(model, res) %>%
    tidyr::unnest(res) %>%
    dplyr::full_join(bio_ref %>% dplyr::select(ID = GO_ID, Biodomain, Subdomain),
                     by = "ID") %>%
    dplyr::full_join(lab_ref, by = c("Biodomain" = "domain")) %>%
    dplyr::filter(!is.na(model)) %>%
    dplyr::mutate(
      Biodomain = dplyr::if_else(is.na(Biodomain), "none", Biodomain),
      Subdomain = dplyr::if_else(is.na(Subdomain), Biodomain, Subdomain)
    ) %>%
    dplyr::mutate(n_sig = length(unique(ID)), .by = Biodomain) %>%
    dplyr::mutate(Biodomain = forcats::fct_reorder(Biodomain, n_sig)) %>%
    dplyr::arrange(Biodomain, p.adjust)
}


#' GSEA BioDomain Distribution Violin + Jitter Plot
#'
#' @param data Output from enrichGSEA_BD_Protein()
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
    ggplot2::scale_x_continuous(limits = c(-ran, ran), breaks = seq(-ran, ran, 1)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(1, 6),
                                   name = expression(paste("-log"[10], "(p-adj)"))) +
    ggplot2::facet_wrap(~label) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray95", color = "gray80"),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      axis.text        = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}


#' GSEA Sub-BioDomain Distribution Plot
#'
#' @param data Output filtered to specific Biodomains
#' @param ran  NES axis limits (default: 3)
plot_gsea_subdomain_distribution <- function(data, ran = 3) {
  ggplot2::ggplot(data, ggplot2::aes(x = NES, y = Subdomain)) +
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
    ggplot2::scale_x_continuous(limits = c(-ran, ran), breaks = seq(-ran, ran, 1)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_continuous(range = c(1, 6),
                                   name = expression(paste("-log"[10], "(p-adj)"))) +
    ggplot2::facet_grid(rows = vars(abbr), cols = vars(model),
                        scales = "free", space = "free", switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray95", color = "gray80"),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      axis.text        = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}


# ==============================================================================
# 11. KEGG PATHWAY ENRICHMENT
# ==============================================================================

#' KEGG Pathway Enrichment via compareCluster — with NULL Guards
#'
#' Runs KEGG ORA for up- and down-regulated proteins per model group.
#' Skips plotting gracefully if no pathways are enriched or gene lists are empty.
#'
#' @param gp1         Character vector of group/model names
#' @param DEG.genes   Data frame: group, Symbol, EntrezGene, log2FoldChange
#' @param n           Number of top pathways to show per group (default: 10)
#' @param kegg_org    KEGG organism code: "mmu"=mouse, "hsa"=human (default: "mmu")
#' @param pval_cutoff p-value cutoff for enrichment (default: 0.1)
#' @param min_genes   Minimum genes per group required to attempt enrichment (default: 3)
kegg.fn <- function(gp1, DEG.genes, n = 10, kegg_org = "mmu",
                    pval_cutoff = 0.1, min_genes = 3) {
  
  # Build gene lists, removing NA Entrez IDs
  dat_up <- purrr::map(gp1, ~ DEG.genes %>%
                         dplyr::filter(group == .x, log2FoldChange > 0) %>%
                         dplyr::pull(EntrezGene) %>%
                         na.omit())
  dat_dn <- purrr::map(gp1, ~ DEG.genes %>%
                         dplyr::filter(group == .x, log2FoldChange < 0) %>%
                         dplyr::pull(EntrezGene) %>%
                         na.omit())
  names(dat_up) <- names(dat_dn) <- gp1
  
  # Drop groups with too few genes for meaningful enrichment
  dat_up <- dat_up[sapply(dat_up, length) >= min_genes]
  dat_dn <- dat_dn[sapply(dat_dn, length) >= min_genes]
  
  # Helper: run compareCluster and plot safely
  run_kegg_plot <- function(gene_list, direction_label) {
    if (length(gene_list) == 0) {
      message("KEGG ", direction_label, ": no groups had >= ", min_genes, " genes. Skipping.")
      return(invisible(NULL))
    }
    
    enr <- tryCatch(
      clusterProfiler::compareCluster(gene_list, fun = "enrichKEGG",
                                      pvalueCutoff = pval_cutoff,
                                      organism = kegg_org),
      error = function(e) {
        message("KEGG ", direction_label, " enrichment error: ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
      message("KEGG ", direction_label, ": no significant pathways found at p < ", pval_cutoff, ".")
      return(invisible(NULL))
    }
    
    print(
      enrichplot::dotplot(enr, showCategory = n) +
        ggplot2::ggtitle(paste("KEGG:", direction_label)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                              face = "bold", size = 11, color = "black")
        )
    )
  }
  
  run_kegg_plot(dat_up, "Upregulated")
  run_kegg_plot(dat_dn, "Downregulated")
}