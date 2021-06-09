Heart\_failure\_IB
================
Komarova Margarita

This is R script to RNA-sequencing analysis of data from article DOI:
10.1161/JAHA.120.017091 using DESeq.

This script has two parts: 1. Comparison healthy donors (HD) and
patients with heart failure (HF) 2. Comparison two groups of patients:
HF\_1 and HF\_2. The second group (HF\_2) is characterized by a lower
(almost critical) fraction compared to the first group (HF\_1).

### Counts table download.

``` r
#setwd("C:/R/IB/Projects/Heart_failure_IB")
data_HF <- read.table('Counts_Rmatrix.txt', sep = '\t', header = T, row.names=1)
```

### Experiment condition table download. Compare healthy donors (HD) and pations with chronic heart failure (HF).

``` r
samples = names(data_HF)
colData = read.table('Condition.txt', sep = '\t', header = T, row.names = 1)
```

### 1\. DESeq2 analysis. Healthy donors (HD) vs pations with heart failure (HF).

``` r
#library(DESeq2)

dds = DESeqDataSetFromMatrix(countData=data_HF, colData=colData, design = ~Condition)
```

Compare patient with healthy donors.

``` r
#Set the reference to be compared
dds$Condition = relevel(dds$Condition,"HD")
```

Run DESeq2.

``` r
dds = DESeq(dds)

# Results
res = results(dds)
res
```

    ## log2 fold change (MLE): Condition HF vs HD 
    ## Wald test p-value: Condition HF vs HD 
    ## DataFrame with 44587 rows and 6 columns
    ##                     baseMean     log2FoldChange             lfcSE
    ##                    <numeric>          <numeric>         <numeric>
    ## DDX11L1                    0                 NA                NA
    ## WASH7P      14.9235647556513  0.271225362581121 0.489338107830215
    ## MIR6859-1   1.48504340514308  0.286350343965533  1.92075189518419
    ## MIR1302-2HG                0                 NA                NA
    ## MIR1302-2                  0                 NA                NA
    ## ...                      ...                ...               ...
    ## TRNS2                      0                 NA                NA
    ## TRNL2       45.8119421954745   0.41909676062714 0.396942540887963
    ## TRNE         5716.7736361353 -0.182162006911235 0.296641372375778
    ## TRNT        167.473790443736  0.194380723793235 0.271139061762957
    ## TRNP        197.207481499932  0.410386869592659 0.318229921729017
    ##                           stat            pvalue              padj
    ##                      <numeric>         <numeric>         <numeric>
    ## DDX11L1                     NA                NA                NA
    ## WASH7P       0.554269856038329 0.579394176830988 0.996788728150696
    ## MIR6859-1    0.149082421672203 0.881488596349595 0.996788728150696
    ## MIR1302-2HG                 NA                NA                NA
    ## MIR1302-2                   NA                NA                NA
    ## ...                        ...               ...               ...
    ## TRNS2                       NA                NA                NA
    ## TRNL2          1.0558121578242 0.291054042575056 0.996788728150696
    ## TRNE        -0.614081594392291 0.539161410188537 0.996788728150696
    ## TRNT         0.716904169135071  0.47343322931475 0.996788728150696
    ## TRNP          1.28959234054087 0.197192236515236 0.996788728150696

Data sorting with p-adj and Log2FoldChange.

``` r
sorted = res[with(res, order(padj, -log2FoldChange)), ]
sorted_df <- data.frame("id"=rownames(sorted),sorted)
sorted_df <- na.omit(sorted_df) 
length(unique(sorted_df$id)) #24631 unique genes
```

    ## [1] 24631

``` r
write.table(sorted_df$id, "sorted.txt", sep = "\t", row.names = F, col.names = F)
```

PCA plot of all samples (Figure 1).

``` r
#PCA plot

png("PCA_plot_HD_HF.png", 500, 500, pointsize=30)
vst_dds <- vst(dds)
counts.norm <- assay(vst_dds)
PCA_plot <- plotPCA(vst_dds,intgroup= c("Condition"))
PCA_plot
dev.off()
```

    ## png 
    ##   2

``` r
PCA_plot
```

![](HF_IB_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Figure 1. PCA plot of sample distribution.

There is a different between samples within the same group.

Volcanolot (Figure 2).

``` r
#Volcano plot
#library(ggplot2)
gdata <- data.frame(
  x=res$log2FoldChange,
  y=-log10(res$pvalue)
)
png("Volcanoplot_HD_HF.png", 500, 500, pointsize=30)
plot <- ggplot(data=gdata, aes(x=x, y=y, color = (res$padj<0.05) )) +
  geom_point(size=1) + theme_bw()  +
  xlab("Log fold change") +
  ylab("Adjusted p.value")
plot
dev.off()
```

    ## png 
    ##   2

``` r
plot
```

![](HF_IB_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Figure 2. Volcanoplot of genes after DESeq2 analysis. Genes with
padj\<0.05 are marked in blue.

There are only two differentially expressed genes: AMPD3 и LDHC.

Next, was done signal pathway analysis using GO analysis.

GO analysis on stat variable.

``` r
#library(data.table)
#library(AnnotationDbi)
sorted_df <- as.data.table(sorted_df)
p_adj_below_0.05_HD_HF <- sorted_df %>% filter(sorted_df$padj < 0.05 & abs(sorted_df$log2FoldChange)>1)
sorted_df$id <- as.character(sorted_df$id)
sorted_df[, row_Genes := mapIds(org.Hs.eg.db, keys=id, 
                                        keytype="SYMBOL", column="ENTREZID")]
# You can saw "keyType" list by "keytypes(org.Hs.eg.db)" command.

sorted_df <- unique(sorted_df, by = "row_Genes")
sorted_df <- na.omit(sorted_df)
original_gene_list1 <- sorted_df$stat
names(original_gene_list1) <- sorted_df$row_Genes
gene_list1 <- na.omit(original_gene_list1)
gene_list1 <- sort(gene_list1, decreasing = TRUE)
gse <- gseGO(geneList=gene_list1, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = FALSE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
```

Number of up и down regulated signal pathways:

``` r
gse_res <- gse@result
gse_down <- gse_res %>% filter(gse_res$NES < 0) #114 down pathways
gse_up <- gse_res %>% filter(gse_res$NES > 0) #156 up pathways
nrow(gse_down)
```

    ## [1] 122

``` r
nrow(gse_up)
```

    ## [1] 153

There are 114 down regulated and 156 up regulated pathways.

GO analysis visualization. Dotplot graph (Figure 3).

``` r
require(DOSE)
png("GO_1.png", 1500, 1000)
GO_plot <- dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
GO_plot
dev.off()
```

    ## png 
    ##   2

``` r
GO_plot
```

![](HF_IB_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Figure 3. Dotplot of Go analysis. HD vs HF.

FGSEA analysis (Figure 4 and 5).

``` r
m_df_stats <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
pathways_stats <- split(m_df_stats$entrez_gene, m_df_stats$gs_name)
fgseaRes_stats <- fgsea(pathways = pathways_stats, 
                  stats    = gene_list1,
                  nperm = 100000, 
                  nproc=4,
                  minSize  = 15,
                  maxSize  = 500)

topPathwaysUp1_stats <- fgseaRes_stats[ES > 0 & padj < 0.05][head(order(padj), n=20), pathway]
topPathwaysDown1_stats <- fgseaRes_stats[ES < 0 & padj < 0.05][head(order(padj),n=20), pathway]
topPathways1_stats <- c(topPathwaysUp1_stats, rev(topPathwaysDown1_stats))
```

``` r
png("FGSEA_1.png", 2500, 1000)
FGSEA_plot <- plotGseaTable(pathways_stats[topPathways1_stats], gene_list1, fgseaRes_stats, gseaParam=0.5)
FGSEA_plot
```

    ## NULL

``` r
dev.off()
```

    ## png 
    ##   2

``` r
FGSEA_plot
```

    ## NULL

Figure 4. FGSEA result visualization. The first 20 up and down regulated
signal pathways.

The other way to FGSEA analysis visualization. The first 20 up and down
regulated signal pathways.

``` r
fgseaRes_up1 <- fgseaRes_stats[ES > 0 & padj < 0.05][head(order(padj,decreasing = FALSE), n=20),]
fgseaRes_down1 <- fgseaRes_stats[ES < 0 & padj < 0.05][head(order(padj,decreasing = FALSE), n=20),]
fgseaRes_up_down1 <- rbind(fgseaRes_up1, fgseaRes_down1)
```

``` r
png("FGSEA_2.png", 1500, 1000)
FGSEA_plot2 <- ggplot(fgseaRes_up_down1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA. HD vs HF") + 
  theme_minimal()
FGSEA_plot2
dev.off()
```

    ## png 
    ##   2

``` r
FGSEA_plot2
```

![](HF_IB_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Figure 5. Hallmark pathways NES from GSEA. HD vs HF.

In both cases, we got the same up and down regulated signal pathways in
patients with heart failure.The pathways that control cellular
respiration and protein localisation were downregulated in patients.
Comparison of healthy donors and patients from our laboratory with fgsea
gave many DEGs and slightly different corresponding signaling pathways.
This can be explained by the not very good choice of the control group
in the article and by another muscle.

# 2\. Compare heart failure patient group.

From PCA plot (Figure 1) we can see heart failure patient division into
two groups.Based on the data of the article, the reason for this
division was identified. Patients were statistically significantly
different in their left ventricular ejection fraction. Further,
signaling pathways were analyzed for these two subgroups. The second
group (HF\_2) of patients has a lower left ventricular ejection
fraction.

``` r
#samples2 <- names(data_HF[,4:9])
colData2 <- read.table('Conditions2.txt', sep = '\t', header = T, row.names = 1)
list <- rownames(colData2)
data_HF2 <- data_HF[,list]
data_HF2 <- data_HF2[, c("X3671309", "X3671311", "X3671312", "X3671310", "X3671313", "X3671314")] #swapped columns 
```

``` r
dds2 = DESeqDataSetFromMatrix(countData=data_HF2, colData=colData2, design = ~Condition)
dds2$Condition = relevel(dds2$Condition,"HF_1")
# Run DESeq
dds2 = DESeq(dds2)

# Format the results.
res2 = results(dds2)
res2
```

    ## log2 fold change (MLE): Condition HF 2 vs HF 1 
    ## Wald test p-value: Condition HF 2 vs HF 1 
    ## DataFrame with 44587 rows and 6 columns
    ##                     baseMean     log2FoldChange             lfcSE
    ##                    <numeric>          <numeric>         <numeric>
    ## DDX11L1                    0                 NA                NA
    ## WASH7P      15.0439517450123 -0.818241835603847 0.578688459190803
    ## MIR6859-1   1.49114484839109  -4.02930818479315  2.06905032099427
    ## MIR1302-2HG                0                 NA                NA
    ## MIR1302-2                  0                 NA                NA
    ## ...                      ...                ...               ...
    ## TRNS2                      0                 NA                NA
    ## TRNL2       47.6200407696696   1.04355811769992 0.375116260751459
    ## TRNE        5204.65138536022  0.458597547511339 0.248320484128539
    ## TRNT        166.246625452665 -0.237222913412382 0.257537181643725
    ## TRNP         204.28294489184 -0.366089921585622  0.37091487150372
    ##                           stat              pvalue               padj
    ##                      <numeric>           <numeric>          <numeric>
    ## DDX11L1                     NA                  NA                 NA
    ## WASH7P       -1.41395913916793   0.157373900140612   0.42544763315671
    ## MIR6859-1    -1.94741913423203  0.0514845059316791                 NA
    ## MIR1302-2HG                 NA                  NA                 NA
    ## MIR1302-2                   NA                  NA                 NA
    ## ...                        ...                 ...                ...
    ## TRNS2                       NA                  NA                 NA
    ## TRNL2         2.78195916009985 0.00540318389104192 0.0508403686009685
    ## TRNE          1.84679709014241  0.0647765537775079  0.259091749360433
    ## TRNT        -0.921121027644679     0.3569872430043  0.650165010673513
    ## TRNP        -0.986991759325968   0.323646678459339  0.617723837864338

``` r
sorted2 <- res2[with(res2, order(padj, -log2FoldChange)), ]
sorted_df2 <- data.frame("id"=rownames(sorted2),sorted2)
sorted_df2 <- unique(sorted_df2)
sorted_df2 <- na.omit(sorted_df2)
p_adj_below_0.05_HF_HF <- sorted_df2 %>% filter(sorted_df2$padj < 0.05 & abs(sorted_df2$log2FoldChange)>1)

length(unique(sorted_df2$id)) #15271 unique genes
```

    ## [1] 15271

PCA plot of two group of patient (Figure 6).

``` r
png("PCA_plot.png", 500, 500, pointsize=30)
vst_dds2 <- vst(dds2)
counts.norm <- assay(vst_dds2)
PCA_plot2 <- plotPCA(vst_dds2,intgroup= c("Condition"))
PCA_plot2
dev.off()
```

    ## png 
    ##   2

``` r
PCA_plot2
```

![](HF_IB_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Figure 6. PCA plot of two group of patient.

Volcano plot of two groups of patient (Figure 7).

``` r
gdata2 <- data.frame(
  x=sorted_df2$log2FoldChange,
  y=-log10(sorted_df2$pvalue)
)
png("Volcanoplot_HF_HF.png", 700, 500, pointsize=30)
volcano_HF <- ggplot(data=gdata2, aes(x=x, y=y, color = (sorted_df2$padj<0.05 & abs(sorted_df2$log2FoldChange)>1) )) +
  geom_point(size=1) + theme_bw()  +
  xlab("Log fold change") +
  ylab("Adjusted p.value")
volcano_HF
dev.off()
```

    ## png 
    ##   2

``` r
volcano_HF
```

![](HF_IB_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Figure 7. Volcanoplot of genes of two groups of patient.

There are a lot of differentially expressed genes. In second group of
patient (with lower left ventricular ejection fraction) are more genes
with negative Log fold change, it means that genes are more expressed in
the first patient group.

Heatmap of differentially expressed genes.

``` r
#library("pheatmap")
select2 <- order(rowMeans(counts(dds2,normalized=TRUE)),
                decreasing=TRUE)[1:50]
ntd2 <- normTransform(dds2)
df2 <- as.data.frame(colData(dds2)[, 2])
rownames(df2) <- colnames(ntd2)

pdf('Heatmap2.pdf')
pheatmap(assay(ntd2)[select2,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2)
invisible(dev.off())
```

GO analysis.

``` r
sorted_df2 <- as.data.table(sorted_df2)
sorted_df2$id <- as.character(sorted_df2$id)
sorted_df2[, row_Genes := mapIds(org.Hs.eg.db, keys=id, 
                                        keytype="SYMBOL", column="ENTREZID")]
original_gene_list2 <- sorted_df2$stat
names(original_gene_list2) <- sorted_df2$row_Genes
gene_list2<-na.omit(original_gene_list2)
gene_list2 = sort(gene_list2, decreasing = TRUE)
gse2 <- gseGO(geneList=gene_list2, 
             ont ="BP", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 50, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = FALSE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
```

``` r
gse2_result <- gse2@result
gse2_down <- gse2_result %>% filter(gse2_result$NES < 0) #31 pathways
nrow(gse2_down)
```

    ## [1] 31

``` r
gse2_up <- gse2_result %>% filter(gse2_result$NES > 0) #927 pathways
nrow(gse2_up)
```

    ## [1] 946

There are 32 down regulated and 975 up regulated pathways.

Dotplot graph (Figure 8).

``` r
#require(DOSE)
png("GO_2.png", 1000, 1000)
GO_plot2 <- dotplot(gse2, showCategory=30, split=".sign") + facet_grid(.~.sign)
GO_plot2
dev.off()
```

    ## png 
    ##   2

``` r
GO_plot2
```

![](HF_IB_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Figure 8. Dotplot graph of GO analysis. HF\_1 group vs HF\_2 group.

FGSEA analysis for two groups of patient.

``` r
#library(fgsea)
#install.packages('msigdbr')
#library(msigdbr)
m_df2 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
pathways2 <- split(m_df2$entrez_gene, m_df2$gs_name)
fgseaRes2 <- fgsea(pathways = pathways2, 
                  stats    = gene_list2,
                  nperm = 100000, 
                  nproc=4,
                  minSize  = 15,
                  maxSize  = 500)
```

``` r
fgseaRes2[order(padj)]
```

    ##                                                 pathway         pval
    ##    1:              GO_ACTIN_CYTOSKELETON_REORGANIZATION 5.095541e-05
    ##    2:                  GO_ACTIN_FILAMENT_BASED_MOVEMENT 6.472492e-05
    ##    3:                GO_ACTIN_MEDIATED_CELL_CONTRACTION 5.460900e-05
    ##    4:              GO_ACTOMYOSIN_STRUCTURE_ORGANIZATION 8.315317e-05
    ##    5:      GO_ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME 1.278772e-05
    ##   ---                                                               
    ## 3587:     GO_INTRASPECIES_INTERACTION_BETWEEN_ORGANISMS 9.950378e-01
    ## 3588:                            GO_NCRNA_TRANSCRIPTION 9.953605e-01
    ## 3589: GO_NEGATIVE_REGULATION_OF_CHROMOSOME_ORGANIZATION 9.955483e-01
    ## 3590:  GO_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS 9.966213e-01
    ## 3591:                       GO_IMMUNOGLOBULIN_SECRETION 9.980660e-01
    ##              padj         ES        NES nMoreExtreme size
    ##    1: 0.002542436  0.3993464  1.8030050            0   92
    ##    2: 0.002542436  0.3677838  1.7578605            0  129
    ##    3: 0.002542436  0.4592951  2.1107061            0  102
    ##    4: 0.002542436  0.3450160  1.7200568            0  172
    ##    5: 0.002542436 -0.5964933 -2.1544319            0   77
    ##   ---                                                    
    ## 3587: 0.996103061  0.1489734  0.5364451        31080   32
    ## 3588: 0.996103061 -0.1580044 -0.5945659        81310  102
    ## 3589: 0.996103061 -0.1571215 -0.6004306        82743  115
    ## 3590: 0.996898895 -0.1581405 -0.6439906        89965  207
    ## 3591: 0.998065969  0.1493183  0.4353767        37671   15
    ##                                  leadingEdge
    ##    1:       8829,301,5911,84168,5156,998,...
    ##    2:       29119,1756,6444,781,4626,857,...
    ##    3:       29119,1756,6444,781,4626,857,...
    ##    4:    1490,9475,55704,4208,11155,8829,...
    ##    5:    6430,10521,10181,6431,4670,5936,...
    ##   ---                                       
    ## 3587:      8195,5108,552,1979,1620,57680,...
    ## 3588:    905,51547,51593,6621,5438,54973,...
    ## 3589: 7517,23133,5394,348654,80169,84464,...
    ## 3590:   5893,83933,55183,9984,1655,23514,...
    ## 3591:              7009,90865,10758,5950,933

``` r
topPathwaysUp2 <- fgseaRes2[ES > 0 & padj < 0.05][head(order(padj), n=20), pathway]
topPathwaysDown2 <- fgseaRes2[ES < 0 & padj < 0.05][head(order(padj), n=20), pathway]
topPathways2 <- c(topPathwaysUp2, rev(topPathwaysDown2))
png("FGSEA1_HF.png", 1500, 1000)
FGSEA1_HF <- plotGseaTable(pathways2[topPathways2], gene_list2, fgseaRes2, 
              gseaParam=0.5)
FGSEA1_HF
```

    ## NULL

``` r
dev.off()
```

    ## png 
    ##   2

``` r
FGSEA1_HF
```

    ## NULL

Figure 9. The first 20 up and down regulated signal pathways. HF1\_1 vs
HF\_2 group.

``` r
fgseaRes_up2 <- fgseaRes2[ES > 0 & padj < 0.05][head(order(padj,decreasing = FALSE), n=20),]
fgseaRes_down2 <- fgseaRes2[ES < 0 & padj < 0.05][head(order(padj,decreasing = FALSE),n=20),]
fgseaRes_up_down2 <- rbind(fgseaRes_up2, fgseaRes_down2)
png("FGSEA_2_HF.png", 1500, 1000)
FGSEA_2_HF_plot <- ggplot(fgseaRes_up_down2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA. HF_1 vs HF_2") + 
  theme_minimal()
FGSEA_2_HF_plot
dev.off()
```

    ## png 
    ##   2

``` r
FGSEA_2_HF_plot
```

![](HF_IB_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

Figure 10. Hallmark pathways NES from GSEA. HF\_1 vs HF\_2.

EnrichGO analysis on differentially expressed genes in two groups of
patient. Now let’s perform EnrichGO analysis for up and down regulated
signal pathways separately.

For DEG with Log2FoldChange\<0:

``` r
down_DEG <- sorted_df2 %>% filter(sorted_df2$padj < 0.05 & abs(sorted_df2$log2FoldChange)>1)
down_DEG <- down_DEG %>% filter(down_DEG$log2FoldChange < 0)
down_DEG <- as.data.table(down_DEG)
nrow(down_DEG)
```

    ## [1] 923

``` r
yy_down_DEG <- enrichGO(down_DEG$row_Genes, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05)


png("EnrichGO_HF_down.png", 500, 500)
enrich_down <- dotplot(yy_down_DEG, showCategory=30)
enrich_down
dev.off()
```

    ## png 
    ##   2

``` r
enrich_down
```

![](HF_IB_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

Figure 11. Dotplot of signal pathways for DEG with Log2FoldChange\<0.
HF\_1 vs HF\_2.

``` r
barplot(yy_down_DEG, showCategory=8)
```

![](HF_IB_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Figure 12. Barplot of signal pathways for DEG with Log2FoldChange\<0.
HF\_1 vs HF\_2.

Graph of signal pathways for DEG with Log2FoldChange\<0. HF\_1 vs HF\_2.

``` r
pdf('Graph_down.pdf')
plotGOgraph(yy_down_DEG) 
```

    ## $dag
    ## A graphNEL graph with directed edges
    ## Number of Nodes = 64 
    ## Number of Edges = 105 
    ## 
    ## $complete.dag
    ## [1] "A graph with 64 nodes."

``` r
invisible(dev.off())
```

EnrichGO for DEG with Log2FoldChange\>0:

``` r
up_DEG <- sorted_df2 %>% filter(sorted_df2$padj < 0.05 & abs(sorted_df2$log2FoldChange)>1)
up_DEG <- up_DEG %>% filter(up_DEG$log2FoldChange > 0)
up_DEG <- as.data.table(up_DEG)
nrow(up_DEG)
```

    ## [1] 139

``` r
yy_up_DEG <- enrichGO(up_DEG$row_Genes, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05)


png("EnrichGO_HF_up.png", 500, 500)
enrich_up <- dotplot(yy_up_DEG, showCategory=10)
enrich_up
dev.off()
```

    ## png 
    ##   2

``` r
enrich_up
```

![](HF_IB_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Figure 13. Dotplot of signal pathways for DEG with Log2FoldChange\>0.
HF\_1 vs HF\_2.

``` r
barplot(yy_up_DEG, showCategory=8)
```

![](HF_IB_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Figure 14. Barplot of signal pathways for DEG with Log2FoldChange\>0.
HF\_1 vs HF\_2.

Graph of signal pathways for DEG with Log2FoldChange\>0. HF\_1 vs HF\_2.

``` r
pdf('Graph_up.pdf')
plotGOgraph(yy_up_DEG) 
```

    ## $dag
    ## A graphNEL graph with directed edges
    ## Number of Nodes = 24 
    ## Number of Edges = 35 
    ## 
    ## $complete.dag
    ## [1] "A graph with 24 nodes."

``` r
invisible(dev.off())
```
