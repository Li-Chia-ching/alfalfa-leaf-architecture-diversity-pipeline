# 🧬 Multidimensional Analysis Pipeline for Alfalfa Leaf Morphogenesis & Plant Architecture

## Overview

This repository provides a **publication-ready, defensive-programming R pipeline (v2.0)** for analyzing compound leaf morphogenesis and plant architecture variation in *Medicago sativa* (alfalfa).

Compared to previous versions, **v2.0 explicitly targets journal reproducibility standards** by integrating **raw data transparency and visualization robustness**.

The pipeline integrates:

* Automated data cleaning and validation
* Composite diversity index calculation (Richness S, Shannon H', Pielou J)
* Correlation analysis (Pearson & Spearman)
* Pattern-level statistical summaries
* Publication-grade visualization with raw data overlay
* **Automatic export of figure-specific Source Data (NEW)**

---

## 🔧 Core Features

### 1. Defensive Programming Framework

* Hard-stop validation of required columns
* Detection of numeric coercion failures
* Explicit error reporting for malformed input

### 2. Dynamic Unit Inference

* Automatically distinguishes:

  * Absolute counts
  * Proportions (auto-converted to counts, N = 15)

### 3. Label Harmonization

* Chinese / mixed labels mapped to standardized English ontology

### 4. Statistical Robustness

* SE-safe computation (no NA propagation in small groups)
* Dual correlation system (Pearson + Spearman)

### 5. Publication-Grade Visualization

* Colorblind-friendly palettes (Viridis, blue–orange)
* Unified ggplot2 styling
* **Raw data overlay (jitter) for transparency (NEW)**

### 6. Reproducible Output System

* Timestamped output directory
* Structured outputs (CSV + PDF)
* **Figure-level Source Data export (NEW)**

---

## 📂 Input Requirements

### Required File

```
RawData_20260501.csv
```

### Required Columns

| Column       | Description                     |
| ------------ | ------------------------------- |
| Family       | Family ID                       |
| ID           | Individual ID                   |
| Plant_Height | Height (cm)                     |
| MF_Total     | Multi-foliate rate              |
| Freq_*       | Leaf pattern frequency or count |

Required frequency columns:

* `Freq_Sym_Term`
* `Freq_Sym_Lat`
* `Freq_Asym_Single`
* `Freq_Asym_Lobed`
* `Freq_Asym_Irreg`

Optional:

* `MF_Pattern`

---

## 🚀 Workflow

1. Environment initialization
2. Data validation (strict schema enforcement)
3. Cleaning + label harmonization
4. Unit inference (count vs proportion)
5. Diversity index computation (S, H', J)
6. Correlation analysis
7. Pattern-level aggregation
8. Visualization + Source Data export

---

## 📊 Output Structure

### Data Tables

| File                               | Description                         |
| ---------------------------------- | ----------------------------------- |
| 01_Cleaned_Data_with_Diversity.csv | Full dataset with diversity indices |
| 02_Architecture_MF_Correlation.csv | Correlation results                 |
| 03_Pattern_Descriptive_Stats.csv   | Pattern-level summary               |

### Figures + Source Data (NEW)

| Figure | Output                       | Source Data         |
| ------ | ---------------------------- | ------------------- |
| Fig1   | Scatter (Architecture vs MF) | SourceData_Fig1.csv |
| Fig2   | Pattern means + raw data     | SourceData_Fig2.csv |
| Fig3   | Shannon violin               | SourceData_Fig3.csv |
| Fig4   | Family bubble plot           | SourceData_Fig4.csv |

---

## 📈 Key Design Philosophy (v7.0)

* **Transparency-first visualization**
  → Every summary plot is paired with raw data

* **Reproducibility compliance**
  → All figures have directly traceable Source Data

* **Reviewer-oriented output**
  → Figures meet typical journal statistical expectations

---

## ⚠️ Data Quality Safeguards

Pipeline will **STOP** if:

* Required columns are missing
* Non-numeric values detected in frequency fields

Pipeline will **WARN** if:

* Mixed unit types detected

---

## 🧪 Dependencies

```r
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)
```
