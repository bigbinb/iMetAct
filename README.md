# iMetAct: An integrated systematic inference of metabolic activity for dissecting tumor metabolic preference and tumor-immune microenvironment
Metabolic enzymes are pivotal catalysts governing biochemical reaction rates, yet systematic profiling of their activities remains technically prohibitive. Traditional experimental approaches (e.g., spectrophotometry, microfluidics) suffer from low throughput, high cost, and incompatibility with complex biological matrices, limiting large-scale metabolic pathway analysis.
Here, we present iMetAct, a framework that integrates metabolic-transcription networks with an information propagation strategy to infer enzyme activity from gene expression data. For more information, Users can also visit our web server **[iMetAct](http://www.imetact.com/)**.
# Installation
```{r}
# install need packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("viper")
install.packages("RANKS")
install.packages("igraph")

if (!requireNamespace("devtools", quietly = TRUE)) {
  # if not have "devtools", install "devtools" first.
  install.packages("devtools")
}
devtools::install_github('xiaolab-xjtu/iMETACT')
```
  # Three Options for Tumor Enzyme Activity Calculation
  ## 1. Pre-computed Network (Recommended for Fast Analysis) 
  The simplest method uses our pre-built enzyme regulon network, covering 33 tumor types from TCGA.

  | Abbr | Name |
  | ---- | ---- |
  | ACC | Adrenocortical carcinoma |
   | BLCA | Bladder Urothelial Carcinoma |
   | BRCA | Breast invasive carcinoma |
   | CESC |Cervical squamous cell carcinoma and endocervical adenocarcinoma |
   | CHOL | Cholangiocarcinoma |
   | COAD | Colon adenocarcinoma |
   | DLBC | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma |
   | ESCA | Esophageal carcinoma |
   | GBM | Glioblastoma multiforme |
   | HNSC | Head and Neck squamous cell carcinoma |
   | KICH | Kidney Chromophobe |
   | KIRC | Kidney renal clear cell carcinoma |
   | KIRP | Kidney renal papillary cell carcinoma |
   | LAML | Acute Myeloid Leukemia |
   | LGG | Brain Lower Grade Glioma |
   | LIHC | Liver hepatocellular carcinoma |
   | LUAD | Lung adenocarcinoma |
   | LUSC | Lung squamous cell carcinoma |
   | MESO | Mesothelioma |
   | OV | Ovarian serous cystadenocarcinoma |
   | PAAD | 	Pancreatic adenocarcinoma |
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
   | UCS | Uterine Carcinosarcoma |
   | UVM | Uveal Melanoma |

## 2. Custom Network Integration
  If you prefer to construct your own metabolic network or identify metabolism-related genes, you can:

 - Use your fully customized regulatory network, or
 - Combine your data with our existing network.
## 3. New Cancer Type Support (No Pre-existing Data Required)
  If your cancer type is not among the 33 provided and you lack resources to build a network de novo, we offer tools to:
  - Construct an enzyme regulon network using our framework, tailored to your cancer type.
# 1. Calculating Cancer Metabolic Enzyme Activity Using Constructed Regulon Networks

**iMetAct** enables rapid metabolic enzyme activity analysis through pre-built regulatory networks derived from TCGA RNA-seq data.

## Key Features:

- **Pre-computed Networks**: Includes ready-to-use metabolic enzyme regulatory networks covering all 33 TCGA cancer types, enabling immediate analysis without additional preprocessing.
  
- **Cancer-Specific Predictions**: Directly computes metabolic enzyme activity profiles for specific tumor types

## Example Implementation:

plaintext
ExpressionMatrix.csv format:
- Rows: Genes
- Columns: Samples
- Content: Tumor gene expression profiles

```
eset <- read.csv('./ExpressionMatrix.csv', row.names = 1)
EA_res <- CalEnzymeActWithTCGAregulon(eset = eset, TCGAtype = 'LIHC') 
```

# 2. Constructing User-Specific Metabolic Regulatory Networks and Identifying Metabolism-Related Genes

iMetAct provides a comprehensive metabolic regulatory network integrating:
- Metabolite-protein/enzyme interactions (MPIs) from KEGG, Reactome, Human-GEM, and BRENDA
- Transcriptional regulatory relationships (TRRs) from TRRUST2
- Protein-protein interactions (PPIs) from STRING

> **Network Features**:
> - Contains over 1 million interaction pairs
> - Directed network format (source in first column, target in second column)
> - Supports user customization: create new networks or extend existing ones

## Step 1: Load Data
Three data components are required:
1. Metabolism-related biological interaction network
2. Metabolite list (optional)
3. Metabolic enzyme list (optional)
```{r}
data(MetNetwork)
data(metabolites)
data(MetabolicEnzymes)
# Import and merge your custom network with:
# MetNetwork <- rbind(MetNetwork, your_network)
```
## Step 2: Identify Metabolism-Related Genes
A restart random walk algorithm simulates metabolic information flow using metabolites as network seeds.

  **Parameters**:

  - filter.pct: Top 20% of genes retained (default=0.2)

  - gamma: Restart probability (default=0.7)

  - tmax/eps: Convergence thresholds
```{r}
MetGenes<- getMetGenes(network = MetNetwork,
                       metabolites = metabolites,
                       filter.pct=0.2,
                       gamma = 0.7,
                       tmax = 1000,
                       eps = 1e-10,
                       norm = TRUE)
```
## Step 3: Infer Metabolic Enzyme Regulatory Network
  1. Filters weakly associated genes from expression matrix

  2. Generates input for network inference
```{r}
# Load and process expression data
expression <- as.matrix(read.csv('./ExpressionMatrix.csv', row.names = 1))
expression_filtered <- expression[rownames(expression) %in% MetGenes, ]

# Prepare ARACNe-AP inputs
write.table(data.frame(regulon = MetabolicEnzymes),
            './path_to_ARACNe/regulon.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(expression_filtered,
            './path_to_ARACNe/Exp_matrix.txt',
            row.names = TRUE, col.names = TRUE)
```
Recommended Tool: ARACNe-AP
**[ARACNe-AP](https://github.com/califano-lab/ARACNe-AP)** for network construction

## Step4: Calculate Enzyme Activity 
Uses VIPER's three-tailed enrichment analysis on inferred networks.
```{r}
regulon <- CreatMetRegulon('ARACNeOutputFile.txt', expression_filtered)
EnzymeActivity <- CalEnzymeAct(expression_filtered, regulon)
```
# 3. Calculating User-Specific Tumor Enzyme Activity

Unlike Method 2, this approach eliminates the need for *de novo* metabolic network construction. You only need to calculate the enzyme regulon network.

## Step 1: Prepare Input Files for Enzyme Regulon Network Calculation

```{r}
# Load and process expression data
expression <- as.matrix(read.csv('./ExpressionMatrix.csv', row.names = 1))

# Load required reference data
data(MetabolicEnzymes)  # Note: Corrected spelling from "MetabolicEnzymes" to "MetabolicEnzymes"
data(Metabolism_Re_genes)

# Filter expression matrix using pre-selected genes
expression_filtered <- expression[rownames(expression) %in% Metabolism_Re_genes, ]

# Generate ARACNe-AP input files
write.table(data.frame(regulon = MetabolicEnzymes),
            './path_to_ARACNe/regulon.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(expression_filtered,
            './path_to_ARACNe/Exp_matrix.txt',
            row.names = TRUE, col.names = TRUE)
```
Recommended Tool: 
**[ARACNe-AP](https://github.com/califano-lab/ARACNe-AP)**
* Output files from this step will be used to construct your enzyme regulatory network.
## Step 2: Calculate Enzyme Activity
The metabolic enzyme activities are inferred using VIPER's three-tailed enrichment analysis.
```{r}
# Create regulatory network object
regulon <- CreatMetRegulon('ARACNeOutputFile.txt', expression_filtered)

# Compute enzyme activities
EnzymeActivity <- CalEnzymeAct(expression_filtered, regulon)
```

# Citations
Wang, Binxian et al. “iMetAct: An integrated systematic inference of metabolic activity for dissecting tumor metabolic preference and tumor-immune microenvironment.” Cell reports vol. 44,3 (2025): 115375.  **[doi:10.1016/j.celrep.2025.115375](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00146-9)**
