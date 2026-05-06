# 🧬 Multidimensional Analysis Pipeline for Alfalfa Leaf Morphogenesis & Plant Architecture

## Overview

This repository provides a **publication-ready, defensive-programming R pipeline (v6.0)** for analyzing compound leaf morphogenesis and plant architecture variation in *Medicago sativa* (alfalfa).

The pipeline integrates:

* Automated data cleaning and validation
* Composite diversity index calculation (Richness S, Shannon H', Pielou J)
* Correlation analysis (Pearson & Spearman)
* Pattern-level statistical summaries with SE-safe computation
* High-resolution, journal-quality visualization (Figure 1–4)

---

## 🔧 Core Features

* **Defensive Programming**

  * Strict column validation
  * Automatic detection of malformed numeric fields
  * NA-safe statistical computation

* **Dynamic Unit Inference**

  * Automatically distinguishes between:

    * Absolute counts
    * Proportions (auto-converted to counts, N = 15)

* **Label Harmonization**

  * Chinese and mixed labels mapped to standardized English terms

* **Robust Visualization System**

  * Publication-grade ggplot2 figures
  * Automatic `cairo_pdf` fallback mechanism

* **Reproducible Output Structure**

  * Timestamped result directory
  * Fully traceable CSV outputs and figures

---

## 📂 Input Requirements

### Required File

```
RawData_20260501.csv
```

### Required Columns

| Column       | Description                              |
| ------------ | ---------------------------------------- |
| Family       | Family ID                                |
| ID           | Individual plant ID                      |
| Plant_Height | Plant height (cm)                        |
| MF_Total     | Multi-foliate rate                       |
| Freq_*       | Frequency or count of leaf pattern types |

Required frequency columns:

* `Freq_Sym_Term`
* `Freq_Sym_Lat`
* `Freq_Asym_Single`
* `Freq_Asym_Lobed`
* `Freq_Asym_Irreg`

Optional:

* `MF_Pattern` (auto-filled if missing)

---

## 🚀 Pipeline Workflow

1. **Environment Initialization**

   * Creates timestamped output directory

2. **Data Validation**

   * Ensures required columns exist
   * Detects invalid numeric conversions

3. **Data Cleaning**

   * Removes invalid rows
   * Standardizes categorical labels

4. **Unit Detection**

   * Identifies whether `Freq_*` are counts or proportions

5. **Diversity Calculation**

   * Computes:

     * Richness (S)
     * Shannon index (H')
     * Pielou’s evenness (J)

6. **Statistical Analysis**

   * Pearson correlation
   * Spearman correlation

7. **Aggregation & Summary**

   * Pattern-level mean, SD, SE
   * Shannon index summary

8. **Visualization Output**

---

## 📊 Output Files

### Data Tables

| File                               | Description                            |
| ---------------------------------- | -------------------------------------- |
| 01_Cleaned_Data_with_Diversity.csv | Cleaned dataset with diversity indices |
| 02_Architecture_MF_Correlation.csv | Correlation results                    |
| 03_Pattern_Descriptive_Stats.csv   | Pattern-level summary statistics       |

### Figures

| Figure | Description                                  |
| ------ | -------------------------------------------- |
| Fig1   | Architecture vs. MF scatter plot             |
| Fig2   | Pattern-wise mean comparison (MF + Height)   |
| Fig3   | Shannon diversity violin plot                |
| Fig4   | Family-level genetic segregation bubble plot |

---

## 📈 Figure Interpretation (Brief)

* **Figure 1**: Quantifies the association between plant height and multi-foliate expression
* **Figure 2**: Compares phenotypic means across dominant leaf patterns
* **Figure 3**: Evaluates intra-plant morphological diversity (Shannon H')
* **Figure 4**: Visualizes family-level genetic segregation structure

---

## ⚠️ Data Quality Safeguards

* Stops execution if:

  * Required columns are missing
  * Non-numeric values detected in frequency columns

* Warns if:

  * Mixed unit types detected (potential contamination)

---

## 🧪 Dependencies

```r
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
```

---

## 🧠 Recommended Use Cases

* Alfalfa population phenotyping
* Morphological diversity quantification
* GWAS phenotype preprocessing
* Reviewer-ready figure generation
