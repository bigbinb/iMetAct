# iMetAct: An integrated systematic inference of metabolic activity for dissecting tumor metabolic preference and tumor-immune microenvironment

Metabolic enzymes are pivotal catalysts governing biochemical reaction rates, yet systematic profiling of their activities remains technically prohibitive. Traditional experimental approaches (e.g., spectrophotometry, microfluidics) suffer from low throughput, high cost, and incompatibility with complex biological matrices, limiting large-scale metabolic pathway analysis.

Here we present **iMetAct**, a framework that integrates metabolic-transcription networks with an information propagation strategy to infer enzyme activity from gene expression data. A web server is also available at **[iMetAct](http://www.imetact.com/)**.

## Publication

Wang, Binxian et al. "iMetAct: An integrated systematic inference of metabolic activity for dissecting tumor metabolic preference and tumor-immune microenvironment." *Cell Reports* 2025, 44(3): 115375. [doi:10.1016/j.celrep.2025.115375](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00146-9)

## Requirements

- R >= 4.1.1
- Packages: `igraph` (>= 1.6.0), `RANKS` (>= 1.0), `viper` (>= 1.32.0)

## Installation

```r
# Install Bioconductor manager (needed for viper)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install dependencies with version constraints
BiocManager::install(c("viper"), version = "3.20")
install.packages(c("igraph", "RANKS"))

# Install iMetAct from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("xiaolab-xjtu/iMetAct")
```

### Reproducible environment

To lock dependency versions for reproducibility, create an installation script:

```r
# Pinned versions known to work with iMetAct 1.0.0
install.packages("igraph", version = "2.3.1")
BiocManager::install("viper", version = "1.46.0")
install.packages("RANKS", version = "1.1")
devtools::install_github("xiaolab-xjtu/iMetAct@v1.0.0")
```

Alternatively, use `renv` to snapshot your environment:

```r
renv::init()
renv::install("igraph@2.3.1")
renv::install("viper@1.46.0")
renv::install("RANKS@1.1")
devtools::install_github("xiaolab-xjtu/iMetAct")
renv::snapshot()
```

---

# Three Options for Tumor Enzyme Activity Calculation

## 1. Pre-computed Network (Recommended for Fast Analysis)

The simplest method uses our pre-built enzyme regulon network, covering 33 tumor types from TCGA.

| Abbr | Name |
| ---- | ---- |
| ACC  | Adrenocortical carcinoma |
| BLCA | Bladder Urothelial Carcinoma |
| BRCA | Breast invasive carcinoma |
| CESC | Cervical squamous cell carcinoma and endocervical adenocarcinoma |
| CHOL | Cholangiocarcinoma |
| COAD | Colon adenocarcinoma |
| DLBC | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma |
| ESCA | Esophageal carcinoma |
| GBM  | Glioblastoma multiforme |
| HNSC | Head and Neck squamous cell carcinoma |
| KICH | Kidney Chromophobe |
| KIRC | Kidney renal clear cell carcinoma |
| KIRP | Kidney renal papillary cell carcinoma |
| LAML | Acute Myeloid Leukemia |
| LGG  | Brain Lower Grade Glioma |
| LIHC | Liver hepatocellular carcinoma |
| LUAD | Lung adenocarcinoma |
| LUSC | Lung squamous cell carcinoma |
| MESO | Mesothelioma |
| OV   | Ovarian serous cystadenocarcinoma |
| PAAD | Pancreatic adenocarcinoma |
| PCPG | Pheochromocytoma and Paraganglioma |
| PRAD | Prostate adenocarcinoma |
| READ | Rectum adenocarcinoma |
| SARC | Sarcoma |
| SKCM | Skin Cutaneous Melanoma |
| STAD | Stomach adenocarcinoma |
| STES | Stomach and Esophageal carcinoma |
| TGCT | Testicular Germ Cell Tumors |
| THCA | Thyroid carcinoma |
| THYM | Thymoma |
| UCEC | Uterine Corpus Endometrial Carcinoma |
| UCS  | Uterine Carcinosarcoma |
| UVM  | Uveal Melanoma |

### Usage

Expression matrix format:
- Rows: Genes
- Columns: Samples
- Content: Tumor gene expression profiles

```r
library(iMetAct)

eset <- read.csv("./ExpressionMatrix.csv", row.names = 1)
EA_res <- CalEnzymeActWithTCGAregulon(eset = eset, TCGAtype = "LIHC")
```

## 2. Custom Network Integration

If you prefer to construct your own metabolic network or identify metabolism-related genes:

- Use your fully customized regulatory network, or
- Combine your data with our existing network.

### Step 1: Load data

Three data components are required:

```r
data(MetNetwork)
data(metabolites)
data(MetabolicEnzymes)

# Optionally merge custom network:
# MetNetwork <- rbind(MetNetwork, your_network)
```

### Step 2: Network propagation

A restart random walk algorithm simulates metabolic information flow using metabolites as network seeds.

| Parameter    | Default | Description |
| ------------ | ------- | ----------- |
| `filter.pct` | 0.2     | Top 20% of genes retained |
| `gamma`      | 0.7     | Restart probability |
| `tmax`       | 1000    | Max iterations |
| `eps`        | 1e-10   | Convergence tolerance |

```r
MetGenes <- getMetGenes(
    network    = MetNetwork,
    metabolites = metabolites,
    filter.pct  = 0.2,
    gamma       = 0.7,
    tmax        = 1000,
    eps         = 1e-10,
    norm        = TRUE
)
```

### Step 3: Metabolic regulon construction

Filter weakly associated genes and prepare ARACNe-AP inputs:

```r
expression <- as.matrix(read.csv("./ExpressionMatrix.csv", row.names = 1))
expression_filtered <- expression[rownames(expression) %in% MetGenes, ]

# Prepare ARACNe-AP input files
write.table(data.frame(regulon = MetabolicEnzymes),
            "./path_to_ARACNe/regulon.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(expression_filtered,
            "./path_to_ARACNe/Exp_matrix.txt",
            row.names = TRUE, col.names = TRUE)
```

Run **[ARACNe-AP](https://github.com/califano-lab/ARACNe-AP)** to infer the regulatory network, then:

```r
# Create regulatory network object
regulon <- CreatMetRegulon(network, expression_filtered)

# Compute enzyme activities
EnzymeActivity <- CalEnzymeAct(regulon, expression_filtered)
```

## 3. New Cancer Type Support (No Pre-existing Data Required)

If your cancer type is not among the 33 pre-computed types, construct a cancer-specific enzyme regulon network using our framework without building the base metabolic network.

### Step 1: Prepare input files

```r
expression <- as.matrix(read.csv("./ExpressionMatrix.csv", row.names = 1))

data(MetabolicEnzymes)
data(Metabolism_Re_genes)

# Filter expression matrix using pre-selected metabolism-related genes
expression_filtered <- expression[rownames(expression) %in% Metabolism_Re_genes, ]

# Generate ARACNe-AP input files
write.table(data.frame(regulon = MetabolicEnzymes),
            "./path_to_ARACNe/regulon.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(expression_filtered,
            "./path_to_ARACNe/Exp_matrix.txt",
            row.names = TRUE, col.names = TRUE)
```

### Step 2: Calculate iMetAct scores

Run **[ARACNe-AP](https://github.com/califano-lab/ARACNe-AP)** first, then:

```r
regulon <- CreatMetRegulon(network, expression_filtered)
EnzymeActivity <- CalEnzymeAct(regulon, expression_filtered)
```
