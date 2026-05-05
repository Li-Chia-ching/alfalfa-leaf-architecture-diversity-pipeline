# ==============================================================================
# 项目：紫花苜蓿复叶形态建成与株形变异多维分析 Pipeline (V1.0)
# 架构：Tidyverse 环境下的自动化数据清洗、复合多样性测算与出版级图表输出
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# --- 1. 自动化输出环境配置 ---
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Architecture_LeafPattern_Results_", current_time)
if(!dir.exists(out_dir)) dir.create(out_dir)
message(">>> 初始化工作空间: ", out_dir)

# --- 2. 数据读取与防御性列名校验 ---
file_path <- "RawData_20260501.csv"
if(!file.exists(file_path)) stop("错误：系统未检测到输入文件: ", file_path)

raw_df <- read_csv(file_path, show_col_types = FALSE)

# 强制防错：核心列名必须存在
required_cols <- c("Family", "ID", "Plant_Height", "MF_Total", 
                   "Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", 
                   "Freq_Asym_Lobed", "Freq_Asym_Irreg")
missing_cols <- setdiff(required_cols, colnames(raw_df))
if(length(missing_cols) > 0) {
  stop("错误：输入文件缺失关键列 -> ", paste(missing_cols, collapse = ", "))
}

# 基础清洗与 NA 填补
clean_df <- raw_df %>%
  mutate(across(c(Plant_Height, MF_Total, starts_with("Freq_")), as.numeric)) %>%
  filter(!is.na(Plant_Height) & !is.na(MF_Total))

if(!"MF_Pattern" %in% colnames(clean_df)) clean_df$MF_Pattern <- "Unclassified"
clean_df$MF_Pattern[is.na(clean_df$MF_Pattern) | clean_df$MF_Pattern == ""] <- "Unclassified"

# --- 3. 高鲁棒性：动态单位复合判定 (频率 vs 频数) ---
freq_cols <- c("Freq_Sym_Term", "Freq_Sym_Lat", "Freq_Asym_Single", "Freq_Asym_Lobed", "Freq_Asym_Irreg")

# 备份原始输入数据，保留证据链
clean_df <- clean_df %>%
  mutate(across(all_of(freq_cols), ~ .x, .names = "{.col}_Raw"))

# 鲁棒判据：检查是否所有非 NA 值都是整数 (允许 1e-6 的浮点误差)
is_all_integers <- function(x) {
  vals <- na.omit(x)
  if(length(vals) == 0) return(TRUE)
  all(abs(vals - round(vals)) < 1e-6)
}

all_int_check <- sapply(clean_df[freq_cols], is_all_integers)

if (all(all_int_check)) {
  message(">>> 智能识别：Freq 列数据皆为整数，确认为绝对频数(Counts)。")
} else {
  max_val <- max(clean_df[freq_cols], na.rm = TRUE)
  if (max_val <= 1.05 && max_val > 0) {
    message(">>> 智能识别：Freq列包含小数且最大值 <= 1，确认为频率(Proportion)。自动乘以鉴定叶片数(N=15)转换为频数。")
    clean_df <- clean_df %>%
      mutate(across(all_of(freq_cols), ~ round(.x * 15)))
  } else {
    warning(">>> 警告：Freq 列既包含非整数，且最大值 > 1.05。请检查原始数据格式是否混淆！目前按原始数据强制运行。")
  }
}

# --- 4. 复合发育多样性指数(S, H, J)计算函数 ---
calc_diversity_indices <- function(row_data) {
  counts <- as.numeric(row_data)
  counts[is.na(counts)] <- 0
  n_multi <- sum(counts)
  
  # 野生型或无明确变异记录的植株，赋予 NA 以避免拉低群体多样性均值
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

# 应用多样性算子并合并
div_results <- bind_rows(apply(clean_df[, freq_cols], 1, calc_diversity_indices))
final_df <- bind_cols(clean_df, div_results)

write_csv(final_df, file.path(out_dir, "01_Cleaned_Data_with_Diversity.csv"))
message(">>> 多维表型矩阵已生成 (01_Cleaned_Data_with_Diversity.csv).")

# --- 5. 统计建模与关联分析 (株形 vs. 多叶表达) ---
cor_p <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "pearson")
cor_s <- cor.test(final_df$Plant_Height, final_df$MF_Total, method = "spearman")

cor_res <- data.frame(
  Variable_X = "Plant_Architecture_Height",
  Variable_Y = "MF_Total",
  Pearson_r = round(cor_p$estimate, 4), P_Value_Pearson = signif(cor_p$p.value, 4),
  Spearman_rho = round(cor_s$estimate, 4), P_Value_Spearman = signif(cor_s$p.value, 4)
)
write_csv(cor_res, file.path(out_dir, "02_Architecture_MF_Correlation.csv"))

# 模式描述性统计 (严谨标注有效计算样本量)
pattern_stats <- final_df %>%
  group_by(MF_Pattern) %>%
  summarise(
    Total_N = n(),
    Mean_MF_Total = round(mean(MF_Total, na.rm=TRUE), 4),
    Mean_Architecture_Height = round(mean(Plant_Height, na.rm=TRUE), 2),
    Valid_N_Shannon = sum(!is.na(Shannon_H)), # 清晰展示多样性计算的底层基数
    Mean_Shannon_H = round(mean(Shannon_H, na.rm=TRUE), 4)
  ) %>% arrange(desc(Total_N))
write_csv(pattern_stats, file.path(out_dir, "03_Pattern_Descriptive_Stats.csv"))

# --- 6. 高级数据可视化 ---

# 图 A: 家系级别的株形与多叶遗传分离气泡图
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

ggsave(file.path(out_dir, "04_Family_Genetics_Plot.pdf"), p_family, width = 9, height = 6)

# 图 B: 复叶形态建成多样性 (Shannon_H) 的发育偏倚解析
# 过滤掉野生型及多样性无意义的数据点
shannon_df <- final_df %>% filter(!is.na(Shannon_H) & MF_Pattern != "无主导型" & MF_Pattern != "Unclassified")

p_shannon <- ggplot(shannon_df, aes(x = reorder(MF_Pattern, Shannon_H, FUN=median, na.rm=TRUE), y = Shannon_H, fill = MF_Pattern)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, color = "darkgray") +
  scale_fill_npg() +
  labs(
    title = "Developmental Canalization breakdown: Shannon_H across patterns",
    subtitle = "Higher index indicates elevated phenotypic plasticity / lower stability",
    x = "Main Compound Leaf Pattern",
    y = "Morphological Diversity Index (Shannon_H)"
  ) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "05_Developmental_Diversity_Plot.pdf"), p_shannon, width = 8, height = 6)

message(">>> 全系统自动化计算完成！所有指标图表已成功导出。")
# ==============================================================================
