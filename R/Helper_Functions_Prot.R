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
# Author: Ravi S Pandey
# ==============================================================================


# ==============================================================================
# 1. SETUP & DEPENDENCIES
# ==============================================================================

needed.packages <- c(
  "ComplexHeatmap", "AnnotationDbi", "circlize", "corrplot", "cowplot",
  "clusterProfiler", "dplyr", "data.table", "EnhancedVolcano",
  "enrichplot", "forcats", "ggplot2", "gridExtra", "ggpubr", "ggnewscale",
  "ggplotify", "ggrepel", "Hmisc", "janitor", "purrr",
  "RColorBrewer", "stringr", "scales", "tidyr", "tibble", "UpSetR",
  "gt", "gprofiler2", "readxl", "readr", "ggh4x", "pcaMethods"
)

# Load libraries silently — missing packages will produce a clear error message
suppressPackageStartupMessages({
  lapply(needed.packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0(
        "Required package '", pkg, "' is not installed.\n",
        "Install it with: install.packages('", pkg, "') ",
        "or BiocManager::install('", pkg, "')"
      ))
    }
    library(pkg, character.only = TRUE)
  })
})


# ==============================================================================
# 2. REFERENCE DATA LOADER
# ==============================================================================
# These functions load optional external reference datasets.
# They are only called if the user enables these analyses in the notebook config.

#' Load BioDomain Reference Data
#' @param biodom_path Path to annotated_biodomains .rds file
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

#' ANOVA + Tukey HSD for Differential Abundance
#'
#' For each protein (columns starting at `protein_col_start`), runs a one-way
#' ANOVA across sample types and extracts Tukey HSD post-hoc results.
#'
#' @param data Data frame with a `SampleType` column and protein abundance columns
#' @param protein_col_start Integer. Column index where protein data begins (default: 6)
#' @return Data frame with F-value, raw p-value, adjusted p-value, and fold-change per protein
anova_DEG_SampleType <- function(data, protein_col_start = 6) {

  # Remove proteins where any group is entirely NA (underpowered)
  cols_to_remove <- data %>%
    dplyr::group_by(SampleType) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~all(is.na(.x)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-SampleType) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~any(.x))) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "col", values_to = "is_na") %>%
    dplyr::filter(is_na) %>%
    dplyr::pull(col)

  df_cleaned <- data %>% dplyr::select(-dplyr::any_of(cols_to_remove))

  # Run ANOVA + Tukey for each protein column
  res_list <- lapply(protein_col_start:ncol(df_cleaned), function(i) {
    tryCatch({
      aov_obj <- aov(df_cleaned[[i]] ~ SampleType, data = df_cleaned)
      a_res   <- anova(aov_obj)
      t_res   <- as.data.frame(TukeyHSD(aov_obj)$SampleType)

      c(Protein = colnames(df_cleaned)[i],
        F_Val = a_res$F[1],
        P_Val = a_res$Pr[1],
        p_adj = t_res$`p adj`,
        diff  = t_res$diff)
    }, error = function(e) NULL)
  })

  final_df <- as.data.frame(do.call(rbind, res_list))
  rownames(final_df) <- final_df$Protein
  return(final_df[, -1])
}


# ==============================================================================
# 4. OVERLAP & UPSET ANALYSIS
# ==============================================================================

#' Extract Up/Down Regulated Protein Lists Per Model
#'
#' Generates named lists of up- and down-regulated proteins for each model,
#' suitable for UpSetR visualization. Assigns `all.up` and `all.down`
#' to the global environment for compatibility with UpSetR.
#'
#' @param gp1 Character vector of model/comparison names
#' @param deg_data Data frame with columns: Symbol, model, diff
overlap.fn <- function(gp1, deg_data) {
  dat.up <- purrr::map(gp1, ~ deg_data %>% dplyr::filter(model == .x, diff > 0) %>% dplyr::pull(Symbol))
  dat.dn <- purrr::map(gp1, ~ deg_data %>% dplyr::filter(model == .x, diff < 0) %>% dplyr::pull(Symbol))
  names(dat.up) <- names(dat.dn) <- gp1
  all.up  <<- dat.up
  all.down <<- dat.dn
}


# ==============================================================================
# 5. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Differential Abundance Barplot
#'
#' @param data Data frame with columns: comparison, DEPs, Direction, Sex
#' @param pval_cutoff Numeric. Used in subtitle label (from global config)
#' @param lfc_cutoff  Numeric. Used in subtitle label (from global config)
barplot_fn <- function(data, pval_cutoff = 0.05, lfc_cutoff = 0) {
  my_colors <- c("Up" = "#D55E00", "Down" = "#0072B2")

  ggplot(data, aes(x = comparison, y = DEPs, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    facet_wrap(~Sex, scales = "free_y") +
    scale_fill_manual(values = my_colors) +
    labs(
      title    = "Differentially Abundant Proteins (DAPs)",
      subtitle = paste0("adjP < ", pval_cutoff, " | |Log2FC| > ", lfc_cutoff),
      x = "", y = "Number of DAPs", fill = "Direction"
    ) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold", size = 11, color = "black"),
      axis.text.y     = element_text(face = "bold", size = 11, color = "black"),
      axis.title      = element_text(face = "bold", size = 13),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle   = element_text(hjust = 0.5, size = 12, face = "italic"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text      = element_text(face = "bold", size = 13),
      legend.position = "top",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}


#' Protein Fold-Change Expression Plot
#'
#' Bar chart of log2 fold-change for selected proteins across comparisons,
#' with significance asterisks for entries below the p-value cutoff.
#'
#' @param deg_df  Data frame of DEG results filtered to one sex
#' @param proteins Character vector of protein symbols to highlight
#' @param sex_label String label used in plot title (e.g., "Females")
#' @param pval_cutoff Numeric. Significance threshold for asterisk display
plot_protein_expression <- function(deg_df, proteins, sex_label, pval_cutoff = 0.05) {

  plot_data <- deg_df %>%
    dplyr::filter(Symbol %in% proteins) %>%
    dplyr::mutate(significant = padj < pval_cutoff) %>%
    dplyr::filter(!is.na(model))

  if (nrow(plot_data) == 0) {
    message("No data found for the requested proteins: ", paste(proteins, collapse = ", "))
    return(invisible(NULL))
  }

  ggplot(plot_data, aes(x = model, y = diff, fill = model)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(significant, "*", "")),
      vjust = ifelse(plot_data$diff > 0, -0.2, 1.2),
      size = 8, fontface = "bold"
    ) +
    facet_wrap(~Symbol, scales = "free_y") +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    labs(
      title    = paste("Protein Abundance Changes:", sex_label),
      subtitle = paste0("* = adjP < ", pval_cutoff),
      x        = "Comparison",
      y        = "log2 Fold Change",
      fill     = "Comparison"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      strip.background = element_rect(fill = "gray90"),
      strip.text      = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.text     = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}


#' Raw Abundance Boxplot for Selected Proteins
#'
#' @param traits_df Data frame with Genotype, Sex, and protein abundance columns
#' @param proteins  Character vector of protein symbols to plot
#' @param group_col Name of the grouping column (default: "Genotype")
plot_protein_raw_abundance <- function(traits_df, proteins, group_col = "Genotype") {

  df_long <- traits_df %>%
    dplyr::select(dplyr::any_of(c("animalName", group_col, "Sex", "Age")),
                  dplyr::all_of(proteins)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(proteins),
                        names_to  = "Protein",
                        values_to = "LogAbundance")

  ggplot(df_long, aes(x = .data[[group_col]], y = LogAbundance,
                      fill = .data[[group_col]])) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 16,
                 outlier.size = 2, notch = FALSE) +
    ylab("log2(abundance)") +
    facet_wrap(~Sex + Protein, scales = "free") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x     = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
      strip.background = element_rect(fill = "gray95"),
      strip.text      = element_text(face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Spectral")
}


# ==============================================================================
# 6. CO-EXPRESSION MODULE CORRELATION ANALYSIS
# ==============================================================================
# These functions correlate mouse model fold-changes with human AD co-expression
# module signatures (e.g., Emory TMT proteomics modules).
# Requires external reference data loaded via load_module_refs().

#' Correlate Mouse DEG Results with Human AD Co-expression Modules
#'
#' @param deg_df     Data frame of DEG results (Symbol, model, diff)
#' @param human_ref  Human reference with columns: Symbol, Module, AD.Control
#' @param sex_label  String for informational messages
#' @return Data frame of Pearson correlations per Module x Model
run_module_correlation <- function(deg_df, human_ref, sex_label = "") {
  deg_df %>%
    dplyr::select(Symbol, model, diff) %>%
    dplyr::inner_join(human_ref, by = "Symbol") %>%
    dplyr::group_by(Module, model) %>%
    dplyr::filter(dplyr::n() >= 3) %>%  # Minimum proteins for valid correlation
    tidyr::nest(data = c(Symbol, diff, AD.Control)) %>%
    dplyr::mutate(
      cor_test    = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["AD.Control"]], method = "pearson")),
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
#' @param cor_data   Output from run_module_correlation()
#' @param ran        Numeric. Color scale limit (e.g., 0.7); auto-computed if NULL
#' @param order_model Character vector defining row order of models
#' @param module_order Character vector defining column order of modules
#' @param title      String. Plot title
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

  ggplot(plot_df, aes(x = Module, y = model)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation))) +
    geom_point(data = dplyr::filter(plot_df, significant),
               color = "black", shape = 0, size = 10, stroke = 1.2) +
    scale_x_discrete(position = "top") +
    scale_size_continuous(range = c(2, 8), guide = "none", limits = c(0, ran)) +
    scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E",
      midpoint = 0, limits = c(-ran, ran),
      name = "Correlation (r)",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15)
    ) +
    labs(x = "Human Brain Proteomics Module", y = NULL, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0, size = 10, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}


#' Scatter Plot of Module Correlations Per Model
#'
#' @param deg_df       DEG results (Symbol, model, diff)
#' @param human_ref    Human reference module data
#' @param mods_to_plot Character vector of module names to visualize
#' @param sex_label    Label for y-axis
plot_module_scatter_correlations <- function(deg_df, human_ref, mods_to_plot, sex_label) {
  for (mod in mods_to_plot) {
    plot_data <- deg_df %>%
      dplyr::select(Symbol, model, diff) %>%
      dplyr::inner_join(human_ref, by = "Symbol") %>%
      dplyr::filter(Module %in% mod)

    if (nrow(plot_data) > 0) {
      p <- ggpubr::ggscatter(
        plot_data,
        x = "AD.Control", y = "diff",
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
    }
  }
}


#' Scatter Plot of Module Correlations from Linear Model Estimates
#'
#' @param deg_df       LM results (Symbol, Variant, Estimate)
#' @param human_ref    Human reference module data
#' @param mods_to_plot Character vector of module names to visualize
#' @param sex_label    Label for y-axis
plot_module_effect_scatter_correlations <- function(deg_df, human_ref, mods_to_plot, sex_label) {
  for (mod in mods_to_plot) {
    plot_data <- deg_df %>%
      dplyr::select(Symbol, Variant, Estimate) %>%
      dplyr::inner_join(human_ref, by = "Symbol") %>%
      dplyr::filter(Module %in% mod)

    if (nrow(plot_data) > 0) {
      p <- ggpubr::ggscatter(
        plot_data,
        x = "AD.Control", y = "Estimate",
        border = TRUE, label = "Symbol",
        add = "reg.line", conf.int = TRUE,
        xlab = paste("Human AD Module:", mod),
        ylab = paste(sex_label, "Mouse models (Estimate)"),
        font.label = 12, repel = FALSE
      ) +
        ggplot2::theme_bw() +
        ggpubr::stat_cor(method = "pearson") +
        ggpubr::grids(linetype = "dashed") +
        ggplot2::facet_wrap(~Variant, scales = "free")
      print(p)
    }
  }
}


# ==============================================================================
# 6.1 LINEAR MODELING UTILITIES
# ==============================================================================

#' Extract Estimates and P-values from LM Summary Objects
#'
#' @param lm_summaries Named list of `summary(lm(...))` objects, one per protein
#' @param coef_names   Character vector of coefficient names to extract
#' @return Long-format data frame with Symbol, Variant, Estimate, P.Value
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


# ==============================================================================
# 6.2 MODULE CORRELATION FROM LINEAR MODEL ESTIMATES
# ==============================================================================

#' Correlate LM Estimates with Human Module Signatures
#'
#' @param lm_results   Output from extract_lm_stats()
#' @param human_ref    Human reference with Symbol, Module, AD.Control
#' @param ordered_variant Character vector defining factor order of Variant
#' @param module_order Character vector defining factor order of Module
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
      cor_test    = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["AD.Control"]], method = "pearson")),
      estimate    = purrr::map_dbl(cor_test, "estimate"),
      p_value     = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      padj        = p.adjust(p_value, method = "fdr"),
      significant = padj < 0.05,
      significant1 = p_value < 0.05
    ) %>%
    dplyr::select(Module, Variant, correlation = estimate, p_value, padj, significant, significant1)

  if (!is.null(module_order))
    cor_results <- cor_results %>%
      dplyr::mutate(Module = factor(Module, levels = module_order))
  if (!is.null(ordered_variant))
    cor_results <- cor_results %>%
      dplyr::mutate(Variant = factor(Variant, levels = ordered_variant),
                    Variant = forcats::fct_rev(Variant))
  cor_results
}


#' Correlation Plot for LM Variant x Module
#'
#' @param data  Output from correlate_lm_with_human_modules()
#' @param ran   Color scale limit (auto-computed if NULL)
protein_variant_corrplot_test <- function(data, ran = NULL) {
  if (is.null(ran))
    ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100

  ggplot(data, aes(x = Module, y = Variant)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation))) +
    geom_point(data = dplyr::filter(data, significant),
               shape = 0, size = 10, stroke = 1.2, color = "black") +
    ggh4x::facet_grid2(rows = vars(group), scales = "free", space = "free", switch = "y",
                strip = ggh4x::strip_themed(
                  background_y = ggh4x::elem_list_rect(fill = "gray95"))) +
    scale_x_discrete(position = "top") +
    scale_size_continuous(guide = "none", limits = c(0, ran), range = c(1, 8)) +
    scale_color_gradient2(limits = c(-ran, ran),
                          low = "#85070C", high = "#164B6E", mid = "white",
                          name = "Correlation") +
    theme_minimal() +
    theme(
      strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),
      axis.text.x       = element_text(angle = 90, hjust = 0, size = 11),
      axis.text.y       = element_text(size = 11),
      panel.grid        = element_blank(),
      legend.position   = "right"
    )
}


# ==============================================================================
# 7. BIODOMAIN ANALYSIS
# ==============================================================================

#' Correlate Model DEG Fold-Changes with Human BioDomain Reference
#'
#' @param deg_fc   DEG results with columns: Symbol, model, diff
#' @param human_ref Human reference with columns: Gene, Term, ampad_fc
#' @param sig      Significance cutoff for p_value
#' @param order_model Character vector for model factor ordering
#' @return Data frame with Pearson correlations per BioDomain x Model
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
      cor_test  = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["ampad_fc"]], method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(biodom = Term, model, correlation = estimate, p_value, significant)

  dat <- dat %>%
    dplyr::mutate(biodom = factor(biodom, levels = domain_order))
  if (!is.null(order_model))
    dat <- dat %>%
      dplyr::mutate(model = factor(model, levels = order_model),
                    model = forcats::fct_rev(model))
  dat
}


#' BioDomain Model Correlation Dot Plot
#'
#' @param data  Output from corr_function_biodomain_protein()
biodom_model_corrplot <- function(data) {
  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100

  ggplot(data, aes(x = biodom, y = model)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation)), alpha = 0.8) +
    geom_point(data = dplyr::filter(data, significant),
               color = "black", shape = 0, size = 10, stroke = 1.1) +
    scale_x_discrete(position = "top") +
    scale_size_continuous(range = c(2, 8), guide = "none", limits = c(0, ran)) +
    scale_color_gradient2(
      low = "#85070C", mid = "white", high = "#164B6E",
      midpoint = 0, limits = c(-ran, ran),
      name = "Correlation",
      guide = guide_colorbar(ticks = FALSE, barheight = 10, frame.colour = "black")
    ) +
    facet_grid(cols = vars(ref), scales = "free", space = "free", switch = "y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      strip.background  = element_rect(fill = "grey95", color = "black"),
      strip.text.x      = element_text(face = "bold", size = 12),
      strip.text.y.left = element_text(face = "bold", size = 12, angle = 0),
      axis.text.x       = element_text(angle = 90, hjust = 0, size = 11, color = "black"),
      axis.text.y       = element_text(size = 11, color = "black"),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "right"
    )
}


#' Correlate Sub-BioDomain with DEG Fold-Changes
#'
#' @param deg_fc   DEG results (Symbol, model, diff)
#' @param human_ref Sub-biodomain reference (Gene, Term, ampad_fc, Biodomain)
#' @param sig      p-value cutoff
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
      cor_test  = purrr::map(data, ~ cor.test(.x[["diff"]], .x[["ampad_fc"]], method = "pearson")),
      estimate  = purrr::map_dbl(cor_test, "estimate"),
      p_value   = purrr::map_dbl(cor_test, "p.value")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significant = p_value < sig) %>%
    dplyr::select(Biodomain, biodom = Term, model, correlation = estimate, p_value, significant)

  dat <- dat %>% dplyr::mutate(Biodomain = factor(Biodomain, levels = domain_order))
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

  ggplot(data, aes(y = biodom, x = model)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation))) +
    geom_point(data = dplyr::filter(data, significant),
               color = "black", shape = 0, size = 9, stroke = 1.1) +
    scale_x_discrete(position = "top") +
    scale_size(guide = "none", limits = c(0, ran), range = c(1, 8)) +
    scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation",
      guide = guide_colorbar(ticks = FALSE, title.position = "top", title.hjust = 0.5)
    ) +
    facet_grid(rows = vars(Biodomain), scales = "free", space = "free", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
      strip.background  = element_rect(fill = "gray96", color = "white"),
      axis.ticks        = element_blank(),
      axis.text.x       = element_text(angle = 90, hjust = 0, size = 11, color = "black"),
      axis.text.y       = element_text(size = 10, color = "black"),
      panel.grid        = element_blank(),
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
      cor_test  = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["ampad_fc"]], method = "pearson")),
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

  ggplot(data, aes(x = biodom, y = Variant)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation))) +
    geom_point(data = dplyr::filter(data, significant),
               aes(colour = correlation), color = "black", shape = 0, size = 9, stroke = 1.1) +
    scale_x_discrete(position = "top") +
    scale_size(guide = "none", limits = c(0, ran), range = c(1, 8)) +
    scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation",
      guide = guide_colorbar(ticks = FALSE)
    ) +
    facet_grid(rows = vars(group), cols = vars(ref),
               scales = "free", space = "free", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      strip.background  = element_rect(fill = "grey95", color = "black"),
      strip.text.x      = element_text(angle = 0, size = 14, face = "bold"),
      strip.text.y.left = element_text(angle = 0, size = 14, face = "bold"),
      axis.ticks        = element_blank(),
      axis.text.x       = element_text(angle = 90, hjust = 0, size = 14, color = "black"),
      axis.text.y       = element_text(size = 14, color = "black"),
      panel.grid        = element_blank(),
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
      cor_test  = purrr::map(data, ~ cor.test(.x[["Estimate"]], .x[["ampad_fc"]], method = "pearson")),
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


#' Sub-BioDomain Portrait Variant Correlation Plot
variant_corrplot_subbiodom.Potrait <- function(data) {
  ran <- ceiling(max(abs(range(data$correlation, na.rm = TRUE))) * 100) / 100

  ggplot(data, aes(y = biodom, x = Variant)) +
    geom_tile(colour = "black", fill = "white") +
    geom_point(aes(colour = correlation, size = abs(correlation))) +
    geom_point(data = dplyr::filter(data, significant),
               aes(colour = correlation), color = "black", shape = 0, size = 9) +
    scale_x_discrete(position = "top") +
    scale_size(guide = "none", limits = c(0, ran)) +
    scale_color_gradient2(
      limits = c(-ran, ran), low = "#85070C", high = "#164B6E", mid = "white",
      name = "Correlation",
      guide = guide_colorbar(ticks = FALSE)
    ) +
    labs(x = NULL, y = NULL) +
    facet_grid(cols = vars(group), rows = vars(Biodomain_label),
               scales = "free", space = "free", switch = "y") +
    theme(
      strip.text.x      = element_text(angle = 0, size = 14),
      strip.text.y.left = element_text(angle = 0, size = 14),
      strip.background.y = element_rect(fill = "grey95"),
      axis.ticks        = element_blank(),
      axis.text.x       = element_text(angle = 90, hjust = 0, size = 14),
      axis.text.y       = element_text(size = 12),
      panel.background  = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "right"
    )
}


# ==============================================================================
# 8. ENRICHMENT UTILITIES — ORA
# ==============================================================================

#' Gene Ontology Over-Representation Analysis with BioDomain Mapping
#'
#' Runs GO enrichment (ORA) for each model x direction combination, then maps
#' results to AD biological domains.
#'
#' @param dat        DEG data frame (Symbol, diff, padj, model)
#' @param bio_ref    BioDomain reference (from load_biodomain_refs()$biodom)
#' @param lab_ref    Domain color labels (from load_biodomain_refs()$dom_lab)
#' @param org_db     OrgDb object for the organism (e.g., org.Mm.eg.db)
#' @param pval_cutoff Significance threshold for gene selection
#' @param universe   Optional background gene set; auto-generated if NULL
#' @return Data frame of enriched GO terms annotated with BioDomains
enrichORA_BD_Protein <- function(dat, bio_ref, lab_ref, org_db,
                                  pval_cutoff = 0.05, universe = NULL) {
  if (is.null(universe)) {
    message("Generating background universe from OrgDb...")
    universe <- clusterProfiler::bitr(
      unique(as.data.frame(AnnotationDbi::select(
        org_db, keys = AnnotationDbi::keys(org_db, "ENTREZID"),
        columns = "SYMBOL", keytype = "ENTREZID"
      ))$SYMBOL),
      fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db
    )$SYMBOL %>% intersect(unique(dat$Symbol))
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
      Biodomain    = tidyr::replace_na(Biodomain, "none"),
      Subdomain    = tidyr::replace_na(Subdomain, "none"),
      signed_logP  = -log10(p.adjust),
      signed_logP  = ifelse(dir == "dn", -1 * signed_logP, signed_logP)
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
  ggplot(enr.ora, aes(x = signed_logP, y = model)) +
    theme_bw() +
    geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dashed", color = "grey40") +
    geom_violin(data = subset(enr.ora, dir == "up"),
                aes(color = color), scale = "width", alpha = 0.1, fill = "grey90") +
    geom_violin(data = subset(enr.ora, dir == "dn"),
                aes(color = color), scale = "width", alpha = 0.1, fill = "grey90") +
    geom_jitter(aes(size = Count, fill = color),
                color = "black", shape = 21, alpha = 0.5, stroke = 0.2,
                width = 0, height = 0.2) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_continuous(range = c(2, 7)) +
    labs(
      title    = "ORA Results Categorized by Biodomains",
      subtitle = "Directional Enrichment Analysis",
      x = "Signed -log10(p-value)", y = "", size = "Gene Count"
    ) +
    facet_wrap(~Biodomain, scales = "free_x") +
    theme(
      axis.text.y     = element_text(size = 11, face = "bold", color = "black"),
      axis.text.x     = element_text(size = 11, color = "black"),
      axis.title      = element_text(size = 13, face = "bold"),
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text      = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
}


#' Sub-BioDomain ORA Distribution Plot
#'
#' @param data  Filtered ORA data with Biodomain, Subdomain, model columns
#' @param title Plot title
plot_subdomain_distribution <- function(data, title = "Sub-BioDomain ORA Distribution") {
  ggplot(data, aes(x = signed_logP, y = Subdomain)) +
    geom_vline(xintercept = 0, lwd = 0.5, lty = 2, color = "gray50") +
    geom_violin(data = dplyr::filter(data, dir == "up"),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_violin(data = dplyr::filter(data, dir == "dn"),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_jitter(aes(size = Count, fill = color),
                color = "grey20", shape = 21, alpha = 0.4,
                width = 0, height = 0.25) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    facet_grid(rows = vars(Biodomain), cols = vars(model),
               scales = "free_y", space = "free_y", switch = "y") +
    theme_minimal() +
    labs(title = title,
         x = expression(paste(signed, " -log"[10], "(p-adj)")), y = NULL) +
    theme(
      strip.background  = element_rect(fill = "gray90", color = "gray70"),
      strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
      strip.text.x      = element_text(size = 12, face = "bold"),
      axis.text.y       = element_text(size = 10),
      axis.text.x       = element_text(size = 11),
      panel.border      = element_rect(color = "gray90", fill = NA),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom"
    )
}


# ==============================================================================
# 9. GSEA WITH BIODOMAIN ANNOTATION
# ==============================================================================

#' Gene Set Enrichment Analysis with BioDomain Mapping
#'
#' @param dat      DEG data frame (Symbol, diff, model)
#' @param bio_ref  BioDomain reference
#' @param lab_ref  Domain color labels
#' @param org_db   OrgDb object
#' @param seed     Integer. Random seed for reproducibility
#' @return Data frame of GSEA results annotated with BioDomains
enrichGSEA_BD_Protein <- function(dat, bio_ref, lab_ref, org_db, seed = 1111) {
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
    dplyr::full_join(bio_ref %>% dplyr::select(ID = GO_ID, Biodomain, Subdomain), by = "ID") %>%
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
  ggplot(data, aes(x = NES, y = model)) +
    geom_vline(xintercept = 0, lwd = 0.5, lty = 2, color = "gray50") +
    geom_violin(data = dplyr::filter(data, NES > 0),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_violin(data = dplyr::filter(data, NES < 0),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_jitter(aes(size = -log10(p.adjust), fill = color),
                color = "grey20", shape = 21, alpha = 0.4,
                width = 0, height = 0.25) +
    scale_x_continuous(limits = c(-ran, ran), breaks = seq(-ran, ran, 1)) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_continuous(range = c(1, 6),
                          name = expression(paste("-log"[10], "(p-adj)"))) +
    facet_wrap(~label) +
    theme_bw() +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme(
      strip.background  = element_rect(fill = "gray95", color = "gray80"),
      strip.text        = element_text(face = "bold", size = 10),
      axis.text         = element_text(size = 9),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom"
    )
}


#' GSEA Sub-BioDomain Distribution Plot
#'
#' @param data Output filtered to specific Biodomains
#' @param ran  NES axis limits (default: 3)
plot_gsea_subdomain_distribution <- function(data, ran = 3) {
  ggplot(data, aes(x = NES, y = Subdomain)) +
    geom_vline(xintercept = 0, lwd = 0.5, lty = 2, color = "gray50") +
    geom_violin(data = dplyr::filter(data, NES > 0),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_violin(data = dplyr::filter(data, NES < 0),
                aes(color = color), scale = "width", alpha = 0.1, fill = NA) +
    geom_jitter(aes(size = -log10(p.adjust), fill = color),
                color = "grey20", shape = 21, alpha = 0.4,
                width = 0, height = 0.25) +
    scale_x_continuous(limits = c(-ran, ran), breaks = seq(-ran, ran, 1)) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_continuous(range = c(1, 6),
                          name = expression(paste("-log"[10], "(p-adj)"))) +
    facet_grid(rows = vars(abbr), cols = vars(model),
               scales = "free", space = "free", switch = "y") +
    theme_bw() +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme(
      strip.background  = element_rect(fill = "gray95", color = "gray80"),
      strip.text        = element_text(face = "bold", size = 10),
      axis.text         = element_text(size = 9),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom"
    )
}


# ==============================================================================
# 10. KEGG PATHWAY ENRICHMENT
# ==============================================================================

#' KEGG Pathway Enrichment via compareCluster
#'
#' Runs KEGG ORA for up- and down-regulated proteins per model group,
#' then generates dot plots. The KEGG organism code is passed as a parameter
#' so this works with mouse, human, rat, or any supported organism.
#'
#' @param gp1          Character vector of group/model names
#' @param DEG.genes    Data frame with columns: group, Symbol, EntrezGene, log2FoldChange
#' @param n            Integer. Number of top pathways to show per group (default: 10)
#' @param kegg_org     KEGG organism code (default: "mmu" for mouse; "hsa" for human)
#' @param pval_cutoff  Numeric. p-value cutoff for enrichment (default: 0.1)
kegg.fn <- function(gp1, DEG.genes, n = 10, kegg_org = "mmu", pval_cutoff = 0.1) {
  dat_up <- purrr::map(gp1, ~ DEG.genes %>%
                         dplyr::filter(group == .x, log2FoldChange > 0) %>%
                         dplyr::pull(EntrezGene))
  dat_dn <- purrr::map(gp1, ~ DEG.genes %>%
                         dplyr::filter(group == .x, log2FoldChange < 0) %>%
                         dplyr::pull(EntrezGene))
  names(dat_up) <- names(dat_dn) <- gp1

  enr_up <- clusterProfiler::compareCluster(
    dat_up, fun = "enrichKEGG", pvalueCutoff = pval_cutoff, organism = kegg_org
  )
  print(
    enrichplot::dotplot(enr_up, showCategory = n) +
      ggplot2::ggtitle("KEGG: Upregulated") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                          face = "bold", size = 11))
  )

  enr_dn <- clusterProfiler::compareCluster(
    dat_dn, fun = "enrichKEGG", pvalueCutoff = pval_cutoff, organism = kegg_org
  )
  print(
    enrichplot::dotplot(enr_dn, showCategory = n) +
      ggplot2::ggtitle("KEGG: Downregulated") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                          face = "bold", size = 11))
  )
}

# End of Helper_Functions_Prot.R
