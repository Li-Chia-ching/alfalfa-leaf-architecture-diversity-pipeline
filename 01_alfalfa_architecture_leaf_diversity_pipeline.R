# ==============================================================================
# Project: Multidimensional Analysis Pipeline for Compound Leaf Morphogenesis and Plant Architecture Variation in Alfalfa (Medicago sativa)
# Version: V1.0 (Defensive Programming & Publication-Ready Edition)
# Framework: Tidyverse-based automated cleaning -> composite diversity quantification -> publication-grade figure output (Fig 1–4)
# Features: Full English label harmonization, dynamic unit detection, SE-safe computation, Cairo_PDF fallback mechanism
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# --- 1. Automated output environment setup ---
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Architecture_LeafPattern_Results_", current_time)
if(!dir.exists(out_dir)) dir.create(out_dir)
message(">>> [Init] Workspace created: ", out_dir)

# Custom safe plotting function (Cairo fallback mechanism)
safe_ggsave <- function(plot, filename, width, height) {
  filepath <- file.path(out_dir, filename)
  tryCatch({
    ggsave(filepath, plot, width = width, height = height, device = cairo_pdf)
  }, error = function(e) {
    warning(sprintf("cairo_pdf is not available in the current environment. Falling back to base pdf device. File: %s", filename))
    ggsave(filepath, plot, width = width, height = height, device = "pdf")
  })
}

# --- 2. Data import and defensive column validation ---
file_path <- "RawData_20260501.csv"
if(!file.exists(file_path)) stop("FATAL ERROR: Input file not found: ", file_path)

raw_df <- read_csv(file_path, show_col_types = FALSE)

# Strict validation: required columns must exist
required_cols <- c("Family", "ID", "Plant_Height", "MF_Total", 
                   "Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", 
                   "Freq_Asym_Lobed", "Freq_Asym_Irreg")
missing_cols <- setdiff(required_cols, colnames(raw_df))
if(length(missing_cols) > 0) stop("FATAL ERROR: Missing required columns -> ", paste(missing_cols, collapse = ", "))

# Backfill missing pattern column if absent
if(!"MF_Pattern" %in% colnames(raw_df)) raw_df$MF_Pattern <- "Unclassified"

# --- 3. Basic cleaning and invalid character detection ---
# Warning interception: detect implicit NA generated during numeric coercion
clean_df <- raw_df %>%
  mutate(across(c(Plant_Height, MF_Total, starts_with("Freq_")), ~ as.numeric(.x))) 

freq_cols <- c("Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", "Freq_Asym_Lobed", "Freq_Asym_Irreg")

# Check for failed coercion
failed_parse <- purrr::map_lgl(freq_cols, ~ any(is.na(clean_df[[.x]]) & !is.na(raw_df[[.x]])))
if(any(failed_parse)) {
  stop("FATAL ERROR: Non-numeric values detected in Freq columns (e.g., commas or text). Please inspect raw data.")
}

# Remove fully invalid rows and standardize pattern labels to English
clean_df <- clean_df %>%
  filter(!is.na(Plant_Height) & !is.na(MF_Total)) %>%
  mutate(MF_Pattern = case_when(
    MF_Pattern %in% c("无主导型", "Wild_Type") ~ "Wild-Type",
    MF_Pattern %in% c("顶小叶基部着生", "Freq_Sym_Term") ~ "SYM-Term",
    MF_Pattern %in% c("侧小叶基部着生", "Freq_Sym_Lat") ~ "SYM-Lat",
    MF_Pattern %in% c("单生附加小叶", "Freq_Asym_Single") ~ "ASYM-Single",
    MF_Pattern %in% c("小叶深裂型", "Freq_Asym_Lobed") ~ "ASYM-Lobed",
    MF_Pattern %in% c("不规则着生", "Freq_Asym_Irreg") ~ "ASYM-Irreg",
    is.na(MF_Pattern) | MF_Pattern == "" | MF_Pattern == "Unclassified" ~ "Unclassified",
    TRUE ~ MF_Pattern
  ))

# --- 4. High-robustness dynamic unit inference (frequency vs. count) ---
clean_df <- clean_df %>% mutate(across(all_of(freq_cols), ~ .x, .names = "{.col}_Raw"))

is_all_integers <- function(x) {
  vals <- na.omit(x)
  if(length(vals) == 0) return(TRUE)
  all(abs(vals - round(vals)) < 1e-6)
}

if (all(sapply(clean_df[freq_cols], is_all_integers))) {
  message(">>> [Data Check] All Freq columns are integers. Interpreted as absolute counts.")
} else {
  max_val <- max(clean_df[freq_cols], na.rm = TRUE)
  if (max_val <= 1.05 && max_val > 0) {
    message(">>> [Data Check] Freq columns contain decimals with max <= 1.05. Interpreted as proportions. Converting to counts (N = 15)...")
    clean_df <- clean_df %>% mutate(across(all_of(freq_cols), ~ round(.x * 15)))
  } else {
    warning(">>> [Warning] Freq columns contain non-integers and max > 1.05. Potential data contamination detected. Proceeding with caution.")
  }
}

# --- 5. Composite developmental diversity indices (S, H, J) ---
calc_diversity_indices <- function(row_data) {
  counts <- as.numeric(row_data)
  counts[is.na(counts)] <- 0
  n_multi <- sum(counts)
  
  if (n_multi == 0) return(data.frame(N_Mutant_Leaves = 0, Richness_S = 0, Shannon_H = NA_real_, Pielou_J = NA_real_))
  
  p <- counts / n_multi
  p <- p[p > 0]
  S <- length(p)
  H <- -sum(p * log(p))
  J <- ifelse(S > 1, H / log(S), NA_real_)
  
  return(data.frame(N_Mutant_Leaves = n_multi, Richness_S = S, Shannon_H = round(H, 4), Pielou_J = round(J, 4)))
}

div_results <- bind_rows(apply(clean_df[, freq_cols], 1, calc_diversity_indices))
final_df <- bind_cols(clean_df, div_results)
write_csv(final_df, file.path(out_dir, "01_Cleaned_Data_with_Diversity.csv"))
message(">>> [Output] Multidimensional phenotype matrix generated.")

# --- 6. Statistical modeling and report generation (SE-safe computation) ---
cor_p <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "pearson")
cor_s <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "spearman", exact = FALSE)

write_csv(data.frame(
  Variable_X = "Plant_Architecture_Height", Variable_Y = "MF_Total",
  Pearson_r = round(cor_p$estimate, 4), P_Value_Pearson = signif(cor_p$p.value, 4),
  Spearman_rho = round(cor_s$estimate, 4), P_Value_Spearman = signif(cor_s$p.value, 4)
), file.path(out_dir, "02_Architecture_MF_Correlation.csv"))

# Pattern-level summary statistics (with SD NA-safe handling)
pattern_plot_df <- final_df %>%
  group_by(MF_Pattern) %>%
  summarise(
    N = n(),
    Mean_MF = mean(MF_Total, na.rm = TRUE),
    SD_MF = sd(MF_Total, na.rm = TRUE),
    SE_MF = ifelse(N > 1 & !is.na(SD_MF), SD_MF / sqrt(N), 0),
    Mean_Architecture = mean(Plant_Height, na.rm = TRUE),
    SD_Arch = sd(Plant_Height, na.rm = TRUE),
    SE_Architecture = ifelse(N > 1 & !is.na(SD_Arch), SD_Arch / sqrt(N), 0),
    Valid_N_Shannon = sum(!is.na(Shannon_H)), 
    Mean_Shannon_H = round(mean(Shannon_H, na.rm=TRUE), 4)
  ) %>% filter(N >= 5)

write_csv(pattern_plot_df, file.path(out_dir, "03_Pattern_Descriptive_Stats.csv"))

# ==============================================================================
# --- 7. Publication-grade visualization (Figure 1–4) ---
# ==============================================================================
message(">>> [Rendering] Generating high-resolution figures...")

# 🔷 Figure 1: Plant architecture vs. multi-foliate rate (individual-level scatter)
p_fig1 <- ggplot(final_df, aes(x = Plant_Height, y = MF_Total)) +
  geom_point(alpha = 0.6, size = 2, color = "#2C7BB6") +
  geom_smooth(method = "lm", se = TRUE, color = "#D7191C", linewidth = 1) +
  labs(
    title = "Figure 1. Association between Plant Architecture and Multi-foliate Expression",
    subtitle = paste0("Pearson r = ", round(cor_p$estimate, 3), " (p = ", signif(cor_p$p.value, 3), "); ",
                      "Spearman ρ = ", round(cor_s$estimate, 3), " (p = ", signif(cor_s$p.value, 3), ")"),
    x = "Plant Architecture (Height component, cm)", y = "Multi-foliate Rate (MF_Total)"
  ) + theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

safe_ggsave(p_fig1, "Fig1_Architecture_vs_MF_scatter.pdf", 7, 5)

# 🔷 Figure 2: Pattern-wise phenotype means (with SE bars)
p_fig2a <- ggplot(pattern_plot_df, aes(x = reorder(MF_Pattern, Mean_MF), y = Mean_MF, fill = MF_Pattern)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = Mean_MF - SE_MF, ymax = Mean_MF + SE_MF), width = 0.2, linewidth = 0.6) +
  scale_fill_npg() +
  labs(title = "Figure 2A. Multi-foliate Expression across Patterns", x = "Dominant Leaf Pattern", y = "Mean MF_Total ± SE") + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_fig2b <- ggplot(pattern_plot_df, aes(x = reorder(MF_Pattern, Mean_Architecture), y = Mean_Architecture, fill = MF_Pattern)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = Mean_Architecture - SE_Architecture, ymax = Mean_Architecture + SE_Architecture), width = 0.2, linewidth = 0.6) +
  scale_fill_npg() +
  labs(title = "Figure 2B. Plant Architecture across Patterns", x = "Dominant Leaf Pattern", y = "Mean Height (cm) ± SE") + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

safe_ggsave(p_fig2a / p_fig2b, "Fig2_Pattern_Mean_Comparison.pdf", 8, 10)

# 🔷 Figure 3: Shannon diversity violin plot
shannon_df <- final_df %>% filter(!is.na(Shannon_H), MF_Pattern != "Unclassified", MF_Pattern != "Wild-Type")

p_fig3 <- ggplot(shannon_df, aes(x = reorder(MF_Pattern, Shannon_H, median, na.rm = TRUE), y = Shannon_H, fill = MF_Pattern)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_npg() +
  labs(
    title = "Figure 3. Intra-plant Morphological Diversity (Shannon H')",
    subtitle = "Higher H' indicates elevated phenotypic plasticity / lower expressivity stability",
    x = "Dominant Leaf Pattern", y = "Morphological Diversity Index (Shannon H')"
  ) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

safe_ggsave(p_fig3, "Fig3_Shannon_Violin.pdf", 8, 6)

# 🔷 Figure 4: Family-level genetic segregation bubble plot
family_stats <- final_df %>%
  group_by(Family) %>%
  summarise(
    Mean_Arch = mean(Plant_Height, na.rm=TRUE),
    Mean_MF = mean(MF_Total, na.rm=TRUE),
    Dominant_Pattern = names(which.max(table(MF_Pattern))),
    Pop_Size = n()
  )
dropped_families <- sum(family_stats$Pop_Size < 3)
message(sprintf(">>> [Filter] %d small families (n < 3) removed.", dropped_families))

family_stats <- family_stats %>% filter(Pop_Size >= 3)

p_fig4 <- ggplot(family_stats, aes(x = Mean_Arch, y = Mean_MF)) +
  geom_point(aes(size = Pop_Size, fill = Dominant_Pattern), alpha = 0.8, shape = 21, color = "black") +
  scale_fill_npg() +
  labs(
    title = "Figure 4. Genetic Segregation of Architecture and Multi-foliate Traits",
    x = "Plant Architecture Index (Height component, cm)", y = "Mean Multi-foliate Rate",
    size = "Family Size", fill = "Dominant Pattern"
  ) + theme_bw() + theme(text = element_text(size = 12))

safe_ggsave(p_fig4, "Fig4_Family_Genetics_Plot.pdf", 9, 6)

message(">>> [Complete] 🎯 Pipeline execution finished. All data and figures successfully archived.")
# ==============================================================================
