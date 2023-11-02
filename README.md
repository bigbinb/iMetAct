# Introduction of iMETACT
Metabolic enzyme activity plays a critical role in regulating the rate of related metabolic reactions.
However, assessing the activities of all relevant enzymes in metabolic pathways currently presents challenges due to the absence of efficient measurement methods. 
To address these challenges, we have developed a computational workflow named iMetAct. You can also visit our web server **[iMetAct](http://www.imetact.com/)**.
# Install
```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  # 如果不存在，则安装 devtools 包
  install.packages("devtools")
}
devtools::install_github('bigbinb/iMETACT')
```
## Quickly starts
iMETACT quickly starts using the metabolic enzyme regulatory network pre-calculated based on TCGA RNA-seq data.
We have finish the TCGA 33 cancers metabolic network, you can calculated the specific tumor metabolic enzyme activity.
```
eset <- read.csv('./ExpressionMatrix.csv',row.names = 1)# read your gene expression matrix
EA_res <- CalEnzymeActWithTCGAregulon(eset = eset,
                            TCGAtype = 'HCC') # cancer types in TCGA
```
## Proceed step by step

---
### Step1: Load data
A total of three data are loaded, including: 
                      1. metabolism-related biological interaction network (including protein interactions, metabolism-related interactions); 
                      2. metabolite list;
                      3. metabolic enzyme list.
```{r}
data(MetNetwork) # Alternatively, your metabolism-related network can be read
data(metabolites) # metabolites
data(MetabolicEnzymes) # Metabolic Enzymes
```
### Step2: Identify metabolism-related genes 
In this step, a restart random walk algorithm is used to simulate the flow of metabolic information on a large-scale metabolism-related biological network with metabolites as network seeds, and then identify metabolism-related genes.
```{r}
MetGenes<- getMetGenes(network = MetNetwork,
                       metabolites = metabolites,
                       filter.pct=0.2,
                       gamma = 0.6,
                       tmax = 1000,
                       eps = 1e-10,
                       norm = TRUE)
```
### Step3:  Infer metabolic enzyme regulatory network
This step first removes genes with weak associations with metabolism from the metabolic expression matrix. 
The resulting metabolism-related expression matrix was then used to infer the metabolic enzyme regulatory network.
We recommend using ARACNe-AP for inferring and constructing your own network.
**[ARACNe-AP](https://github.com/califano-lab/ARACNe-AP)**

```{r}
# read gene expression matrix
expression <- read.csv('./ExpressionMatrix.csv',row.names = 1)
expression <- as.matrix(expression)
# filter gene expression matrix
expression_filtered <- expression[rownames(expression)%in%MetGenes,]
# write input files for ARACNe-AP
write.table(x = data.frame(regulon=MetabolicEnzymes),
            './path to ARACNe/regulon.txt',
            row.names = F,
            col.names = F)
write.table(x = expression_filtered,
            './path to ARACNe/Exp_matrix.txt',
            row.names = T,
            col.names = T)
# Run ARACNe-AP to construct GRNs


```
### Step4: Calculate enzyme activity 
Finally, we inferred metabolic enzyme activities through viper's three-tailed enrichment analysis method.
```{r}
# creat  metabolic enzyme regulatory network object from ARACNe-AP results
regulon <- CreatMetRegulon('ARACNeOutputFile.txt',
                          expression_filtered)
# infer enzyme activity
EnzymeActivity <- CalEnzymeAct(expression_filtered,
                               regulon)
```
# Citations
xxxxxx
