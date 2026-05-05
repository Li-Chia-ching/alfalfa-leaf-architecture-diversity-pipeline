# ==============================================================================
# Project: Multidimensional Analysis Pipeline for Alfalfa Compound Leaf Morphogenesis 
#          and Plant Architecture Variation (Publication-ready V4.0)
# Framework: Automated data cleaning, composite diversity quantification, and 
#            publication-grade visualization under the Tidyverse ecosystem
# Updates:
#   - Resolved Spearman tie warnings in large datasets
#   - Enabled cairo_pdf device to fix Unicode/Chinese rendering issues in PDF output
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
message(">>> Workspace initialized: ", out_dir)

# --- 2. Data import and defensive column validation ---
file_path <- "RawData_20260501.csv"
if(!file.exists(file_path)) stop("ERROR: Input file not found: ", file_path)

raw_df <- read_csv(file_path, show_col_types = FALSE)

# Strict validation: required core columns must exist
required_cols <- c("Family", "ID", "Plant_Height", "MF_Total", 
                   "Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", 
                   "Freq_Asym_Lobed", "Freq_Asym_Irreg")
missing_cols <- setdiff(required_cols, colnames(raw_df))
if(length(missing_cols) > 0) {
  stop("ERROR: Missing required columns -> ", paste(missing_cols, collapse = ", "))
}

# Basic cleaning and NA handling
clean_df <- raw_df %>%
  mutate(across(c(Plant_Height, MF_Total, starts_with("Freq_")), as.numeric)) %>%
  filter(!is.na(Plant_Height) & !is.na(MF_Total))
# New: Standardize Chinese into Academic English
mutate(MF_Pattern = case_when(
  MF_Pattern == "无主导型" ~ "Wild_Type",
  is.na(MF_Pattern) | MF_Pattern == "" ~ "Unclassified",
  TRUE ~ MF_Pattern
))

if(!"MF_Pattern" %in% colnames(clean_df)) clean_df$MF_Pattern <- "Unclassified"
clean_df$MF_Pattern[is.na(clean_df$MF_Pattern) | clean_df$MF_Pattern == ""] <- "Unclassified"

# --- 3. Robust dynamic unit detection (proportion vs. count) ---
freq_cols <- c("Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", "Freq_Asym_Lobed", "Freq_Asym_Irreg")

# Preserve raw inputs for traceability
clean_df <- clean_df %>%
  mutate(across(all_of(freq_cols), ~ .x, .names = "{.col}_Raw"))

# Robust criterion: check whether all non-NA values are integers (tolerance = 1e-6)
is_all_integers <- function(x) {
  vals <- na.omit(x)
  if(length(vals) == 0) return(TRUE)
  all(abs(vals - round(vals)) < 1e-6)
}

all_int_check <- sapply(clean_df[freq_cols], is_all_integers)

if (all(all_int_check)) {
  message(">>> Auto-detected: Freq columns are integers → treated as absolute counts.")
} else {
  max_val <- max(clean_df[freq_cols], na.rm = TRUE)
  if (max_val <= 1.05 && max_val > 0) {
    message(">>> Auto-detected: Freq columns are proportions (max ≤ 1). Converting to counts using N = 15 leaves.")
    clean_df <- clean_df %>%
      mutate(across(all_of(freq_cols), ~ round(.x * 15)))
  } else {
    warning(">>> WARNING: Mixed or inconsistent Freq format detected (non-integer values with max > 1.05). Please verify input data. Proceeding without transformation.")
  }
}

# --- 4. Composite developmental diversity indices (S, H, J) ---
calc_diversity_indices <- function(row_data) {
  counts <- as.numeric(row_data)
  counts[is.na(counts)] <- 0
  n_multi <- sum(counts)
  
  # Assign NA to wild-type or non-mutant individuals to avoid biasing population-level estimates
  if (n_multi == 0) {
    return(data.frame(N_Mutant_Leaves = 0, Richness_S = 0, Shannon_H = NA_real_, Pielou_J = NA_real_))
  }
  
  p <- counts / n_multi
  p <- p[p > 0]
  S <- length(p)
  H <- -sum(p * log(p))
  J <- ifelse(S > 1, H / log(S), NA_real_)
  
  return(data.frame(
    N_Mutant_Leaves = n_multi,
    Richness_S = S, 
    Shannon_H = round(H, 4), 
    Pielou_J = round(J, 4)
  ))
}

# Apply diversity computation and merge results
div_results <- bind_rows(apply(clean_df[, freq_cols], 1, calc_diversity_indices))
final_df <- bind_cols(clean_df, div_results)

write_csv(final_df, file.path(out_dir, "01_Cleaned_Data_with_Diversity.csv"))
message(">>> Multidimensional phenotypic matrix generated (01_Cleaned_Data_with_Diversity.csv).")

# --- 5. Statistical modeling and association analysis (architecture vs. multi-foliate traits) ---
cor_p <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "pearson")

# Improvement: set exact = FALSE to suppress tie warnings in large samples
cor_s <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "spearman", exact = FALSE)

cor_res <- data.frame(
  Variable_X = "Plant_Architecture_Height",
  Variable_Y = "MF_Total",
  Pearson_r = round(cor_p$estimate, 4), P_Value_Pearson = signif(cor_p$p.value, 4),
  Spearman_rho = round(cor_s$estimate, 4), P_Value_Spearman = signif(cor_s$p.value, 4)
)
write_csv(cor_res, file.path(out_dir, "02_Architecture_MF_Correlation.csv"))

# Pattern-level descriptive statistics (explicit reporting of valid sample sizes)
pattern_stats <- final_df %>%
  group_by(MF_Pattern) %>%
  summarise(
    Total_N = n(),
    Mean_MF_Total = round(mean(MF_Total, na.rm=TRUE), 4),
    Mean_Architecture_Height = round(mean(Plant_Height, na.rm=TRUE), 2),
    Valid_N_Shannon = sum(!is.na(Shannon_H)), 
    Mean_Shannon_H = round(mean(Shannon_H, na.rm=TRUE), 4)
  ) %>% arrange(desc(Total_N))
write_csv(pattern_stats, file.path(out_dir, "03_Pattern_Descriptive_Stats.csv"))

# --- 6. Advanced data visualization ---

# Figure A: Family-level genetic segregation (bubble plot of architecture vs. multi-foliate traits)
family_stats <- final_df %>%
  group_by(Family) %>%
  summarise(
    Mean_Arch = mean(Plant_Height, na.rm=TRUE),
    Mean_MF = mean(MF_Total, na.rm=TRUE),
    Dominant_Pattern = names(which.max(table(MF_Pattern))),
    Pop_Size = n()
  ) %>% filter(Pop_Size >= 3)
write_csv(family_stats, file.path(out_dir, "04_Family_Genetics_PlotData.csv"))

p_family <- ggplot(family_stats, aes(x = Mean_Arch, y = Mean_MF)) +
  geom_point(aes(size = Pop_Size, fill = Dominant_Pattern), alpha = 0.8, shape = 21, color = "black") +
  scale_fill_npg() +
  labs(
    title = "Genetic Segregation of Plant Architecture and Multi-foliate Traits",
    x = "Plant Architecture Index (Height component, cm)",
    y = "Mean Multi-foliate Rate",
    size = "Family Size", fill = "Dominant Pattern"
  ) + theme_bw() + theme(text = element_text(size = 12))

# Improvement: use cairo_pdf to fully resolve Unicode rendering issues
ggsave(file.path(out_dir, "04_Family_Genetics_Plot.pdf"), p_family, 
       width = 9, height = 6, device = cairo_pdf)

# Figure B: Developmental bias in compound leaf morphogenesis (Shannon diversity)
# Exclude wild-type or non-informative samples
shannon_df <- final_df %>% 
  filter(!is.na(Shannon_H) & MF_Pattern != "无主导型" & MF_Pattern != "Unclassified")

p_shannon <- ggplot(shannon_df, aes(x = reorder(MF_Pattern, Shannon_H, FUN=median, na.rm=TRUE), y = Shannon_H, fill = MF_Pattern)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, color = "darkgray") +
  scale_fill_npg() +
  labs(
    title = "Developmental Canalization Breakdown: Shannon_H Across Patterns",
    subtitle = "Higher values indicate increased phenotypic plasticity and reduced developmental stability",
    x = "Dominant Compound Leaf Pattern",
    y = "Morphological Diversity Index (Shannon_H)"
  ) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Improvement: use cairo_pdf to fully resolve Unicode rendering issues
ggsave(file.path(out_dir, "05_Developmental_Diversity_Plot.pdf"), p_shannon, 
       width = 8, height = 6, device = cairo_pdf)

message(">>> Pipeline execution complete. All outputs exported successfully with full Unicode support.")
# ==============================================================================
