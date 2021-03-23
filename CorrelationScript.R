################################################################################
#                                                                              #
#                             TCGA Data Analysis                               #
#                                                                              #
################################################################################

### ============================================================================
### Structure of the script
### ----------------------------------------------------------------------------
###
### 1) TCGA data preparation
### 2) GO analysis of co-correlation genes with Tak1l
### 3) PGE2 production estimation models
### 4) Cox regression analysis
###



### ============================================================================
###
### 1) TCGA data preparation
###
### ============================================================================

### ----------------------------------------------------------------------------
### Imported libraries
### ----------------------------------------------------------------------------
library(RTCGA)
library(RTCGA.clinical)
library(survival)
library(survminer)
library(ggplot2)
library(ggbeeswarm)
library(viridis)
library(gridExtra)

### ----------------------------------------------------------------------------
### Import TCGA expression data
### ----------------------------------------------------------------------------
### *NOTE*
### Data files were downloaded with RTCGA and stored as flat files prior to this
### analysis
expressionDataPath = "./data_RNA_Seq_v2_expression_median.txt"
clinicalDataPath = "./data_bcr_clinical_data_patient.txt"

tcga = read.table(expressionDataPath, header = T)
tcgaExpression = tcga
rownames(tcgaExpression) = tcgaExpression$Entrez_Gene_Id
tcgaExpression = tcgaExpression[, 3:ncol(tcgaExpression)]
colnames(tcgaExpression) = paste0(
    sapply(strsplit(colnames(tcgaExpression), "\\."), `[`, 1),
    "-",
    sapply(strsplit(colnames(tcgaExpression), "\\."), `[`, 2),
    "-",
    sapply(strsplit(colnames(tcgaExpression), "\\."), `[`, 3)
)
# remove paitients with multiple meassurements
tcgaExpression = tcgaExpression[!duplicated(colnames(tcgaExpression))]
# match gene names and entrez ids
tcgaGeneLink = tcga[, c(1, 2)]
# load patient info
tcgaPatientInfo = read.table(
    clinicalDataPath,
    header = T,
    skip = 4,
    fill = T,
    sep = "\t"
)
tcgaPatientInfo$PATIENT_ID = as.character(tcgaPatientInfo$PATIENT_ID)
tcgaExpressionActive = tcgaExpression[rowMeans(tcgaExpression) > 1, ]



### ============================================================================
###
### 2) GO analysis of co-correlation genes with TAK1L
###
### ============================================================================

### ----------------------------------------------------------------------------
### Calculate all pairwise correlations to TAK1L
### ----------------------------------------------------------------------------

# select TAK1L properties
tak1l = subset(tcgaGeneLink, Hugo_Symbol == "MAP3K7CL") # TAK1L
tak1lExpression = log2(tcgaExpression[rownames(tcgaExpression) == tak1l$Entrez_Gene_Id, ] + 0.1) # add pseudocount

# calculate pairwise correlation of TAK1L vs. all expressed genes
corr = apply(tcgaExpressionActive, 1, function(x) {
    corcoef = as.numeric(cor.test(t(tak1lExpression), log2(x + 1),
                                  method = "pearson")$estimate)
    corpval = as.numeric(cor.test(t(tak1lExpression), log2(x + 1),
                                  method = "pearson")$p.value)
    res = cbind.data.frame(corcoef, corpval)
    return(res)
})
corrRes = do.call("rbind", corr)

# format results
corrRes$padj = p.adjust(corrRes$corpval, method = "BH") # multiple testing correction
corrRes$sig = ifelse(corrRes$padj < 0.05, T, F)
corrRes$baseMean = rowMeans(tcgaExpressionActive)
idx = match(rownames(corrRes), tcgaGeneLink$Entrez_Gene_Id)
corrRes$names = tcgaGeneLink$Hugo_Symbol[idx]
corrRes$sigAll = ifelse(
    corrRes$padj < 0.05 & abs(corrRes$corcoef) > 0.3,
    "Strong",
    ifelse(corrRes$padj < 0.05 &
               abs(corrRes$corcoef) <= 0.3, "Weak", "Not")
)

# plotting results
pl = plot.correlationsTak1l(corrRes)
grid.arrange(pl$p1, pl$p2, pl$p3, pl$p4, ncol = 2, nrow = 2)


### ----------------------------------------------------------------------------
### Enrichment analysis of co-correlated subgroup
### ----------------------------------------------------------------------------
###
### calculate all-vs-all correlation among genes highly correlated to TAK1L
corrAll = subset(t(tcgaExpressionActive), select = rownames(corrSig))
corrAll = log2(corrAll + 0.1) # transform to log2 space and add pseudocount
corrAll = cor(corrAll, method = "pearson")

# subset all correlations
goodCorr = corrAll[colnames(corrAll) == tak1l$Entrez_Gene_Id, ]
goodCorr = names(goodCorr[abs(goodCorr) > 0.3])
goodCorr = subset(corrAll, select = goodCorr)
goodCorr = goodCorr[rownames(goodCorr) %in% colnames(goodCorr), ]

# format names
idx = match(rownames(goodCorr), tcgaGeneLink$Entrez_Gene_Id)
rownames(goodCorr) = tcgaGeneLink$Hugo_Symbol[idx]
colnames(goodCorr) = NULL
geneNames = rownames(goodCorr)

# prepare data for GO testing
testDf = data.frame(id = tcgaGeneLink$Entrez_Gene_Id[tcgaGeneLink$Hugo_Symbol %in% geneNames],
                    type = "selected")
background = rownames(corrRes) # make background universe (all expressed genes)

# test enrichment
enrichmentReactome = compareCluster(
    id ~ type,
    data = testDf,
    fun = "enrichPathway",
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    universe = background
)

enrichmentGo = compareCluster(
    id ~ type,
    data = testDf,
    OrgDb = org.Hs.eg.db,
    fun = "enrichGO",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.05,
    universe = background
)

# plotting results
p1 = clusterProfiler::dotplot(
    enrichmentReactome,
    showCategory = 15,
    title = "Reactome pathways",
    font.size = 8
)
p3 = clusterProfiler::dotplot(
    enrichmentGo,
    showCategory = 15,
    title = "GO biological process",
    font.size = 8
)
grid.arrange(p1, p3, ncol = 1)



### ============================================================================
###
### 3) PGE2 production estimation models
###
### ============================================================================

### ----------------------------------------------------------------------------
### PDGFRB-based regression model to capture PGE2 production
### ----------------------------------------------------------------------------

# relevant gene names
nameProd = c("PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3", "HPGD")
nameFibro = c("FAP", "COL1A1", "ACTA2", "PDGFRB", "S100A4", "SERPINH1")

# construct dataframe for PGE2-producing genes
link = tcgaGeneLink[tcgaGeneLink$Hugo_Symbol %in% c(nameProd, "PDGFRB"), ]
dataProd = tcgaExpressionActive[rownames(tcgaExpressionActive) %in% link$Entrez_Gene_Id, ] %>% t
idx = match(colnames(dataProd), link$Entrez_Gene_Id)
colnames(dataProd) = link$Hugo_Symbol[idx]

# compute regression model
dataProd = as.data.frame(log2(dataProd + 0.1)) # add pseudocount
modProd = lm(formula = PDGFRB ~ PTGS1 + PTGS2 + HPGD + PTGES + PTGES2 + PTGES3,
             data = dataProd)
combRes = cbind(dataProd, mPDGFRB = (modProd$fitted.values) * -1) # muliply with -1 to make the model caputre production, rather then degradation
combRes = combRes[, -c(2)]

# reorder
colOrder = c("mPDGFRB", nameProd)
combRes = combRes[, colOrder]

# make plot
ggpairs(
    combRes,
    upper = list(continuous = pairsCorColor),
    lower = list(continuous = pairsBsScatterRast),
    diag = list(continuous = pairsDensity),
    progress = FALSE
) + ggtitle("PDGFRB model")


### ----------------------------------------------------------------------------
### PDGFRB-based regression model agrees with other fibroblast markers
### ----------------------------------------------------------------------------

# construct dataframe for fibroblast markers
link = tcgaGeneLink[tcgaGeneLink$Hugo_Symbol %in% c(nameFibro), ]
dataFibro = tcgaExpressionActive[rownames(tcgaExpressionActive) %in% link$Entrez_Gene_Id, ] %>% t
idx = match(colnames(dataFibro), link$Entrez_Gene_Id)
colnames(dataFibro) = link$Hugo_Symbol[idx]

dataFibro = as.data.frame(log2(dataFibro + 0.1))
combFibro = cbind(dataFibro, mPDGFRB = (modProd$fitted.values) * -1) # muliply with -1 to make the model caputre production, rather then degradation

# reorder
colOrder = c("mPDGFRB", nameFibro)
combFibro = combFibro[, colOrder]

ggpairs(
    combFibro,
    upper = list(continuous = pairsCorColor),
    lower = list(continuous = pairsBsScatterRast),
    diag = list(continuous = pairsDensity),
    progress = FALSE
) + ggtitle("Fibroblas marker")



### ============================================================================
###
### 4) Cox regression analysis
###
### ============================================================================

### ----------------------------------------------------------------------------
### Prepare clinical dataset
### ----------------------------------------------------------------------------

# download total clinical data for breast cancer dataset
clin = survivalTCGA(
    BRCA.clinical,
    extract.cols = c(
        "patient.gender",
        "patient.days_to_birth",
        "patient.vital_status"
    )
)
# match total data with present data
idx = match(colnames(tcgaExpression), clin$bcr_patient_barcode)
clinData = clin[idx, ]
colnames(clinData) = c("times", "id", "status", "gender", "age")
clinData$age = as.numeric(clinData$age) * -1
selProd = tcgaGeneLink[tcgaGeneLink$Hugo_Symbol %in% nameProd, ]

targetExpression = t(as.data.frame(tcgaExpressionActive[
    rownames(tcgaExpressionActive) %in% selProd$Entrez_Gene_Id, ]))
idx = match(colnames(targetExpression), selProd$Entrez_Gene_Id)
colnames(targetExpression) = selProd$Hugo_Symbol[idx]
targetExpression = targetExpression[, nameProd]

# plot expression levels
df = melt(targetExpression)
df$class = ifelse(
    df$Var2 == "HPGD",
    "Degradation",
    ifelse(
        df$Var2 == "PTGS2" |
            df$Var2 == "PTGS1",
        "Production 1",
        "Production 2"
    )
)

p1 = ggplot(df, aes(
    x = Var2,
    y = log2(value),
    color = class
)) +
    geom_beeswarm_rast(
        size = 0.8,
        priority = "density",
        cex = .35,
        alpha = 0.5
    ) +
    theme_bw() +
    geom_boxplot(
        outlier.color = NA,
        color = "black",
        alpha = 0,
        size = 0.5
    ) +
    # scale_color_brewer(palette = "Dark2") +
    scale_color_npg() +
    theme(legend.position = "right") +
    labs(tag = "A") +
    xlab(NULL) + ylab("Expression level [log2(TPM)]")

# calculate and plot Cox model
testData = cbind(clinData, log2(targetExpression + 0.01))
mod = coxph(Surv(times, status) ~ PTGES + PTGES2 + PTGES3 + PTGS1 + PTGS2 + HPGD,
            data = testData)

p2 = ggforest(mod, data = testData, main = "Cox Model (HR)") +
    labs(tag = "B")

grid.arrange(p1, p2, heights = c(1, 1))
