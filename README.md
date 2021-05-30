# Transcriptomic_analysis_of_patients_with_heart_failure
This is spring project in Bioinformatic Institute 2021

Author: Komarova Margarita

In heart failure (HF), functional and metabolic alterations are detected not only in cardiac muscle but also in skeletal muscle tissue. In our laboratory based on Almazov National Medical Research Centre we have done transcriptomic analysis of calf muscle of healthy donors and patients with heart failure. Then we have found the other article (DOI: 10.1161/JAHA.120.017091) where transcriptome analysis of the pectoral muscle was also performed. 
**Aim**
The aim of this project was to investigate the muscle transcriptome of patients with heart failure in an article and compare these results with the same experiment performed in our laboratory.

## Methods:
- FastQC (v0.11.9)
- STAR (v2.7.3a) alignment
- featureCounts (v2.0.1) (getting counts table)
- DESeq2 (v1.26.0) (gene expression analysis in Rstudio)
- also we use fgsea (v1.12.0) and ClusterProfiler (v3.14.3) packages to identify signal pathways. 

This project has two parts:
1. Comparison healthy donors (HD) and patients with heart failure (HF)
2. Comparison two groups of patients: HF_1 and HF_2. The second group (HF_2) is characterized by a lower (almost critical) left ventricular ejection fraction compared to the first group (HF_1).

## The first step of all piplines was quality assessment of reads using FastQC program (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Reads of this article were posted in open access and are available via the link https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8531/samples/. The reads have of high quality and do not need trimming procedure.

## At the next step, the genome was indexed in the STAR program. Example of command for this step:
```
/home/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /mnt4/transciptome_compare/komarova/STAR_genome_index --genomeFastaFiles /mnt4/transciptome_compare/komarova/GCF_000001405.39_GRCh38.p13_genomic.fna --sjdbGTFfile /mnt4/transciptome_compare/komarova/GCF_000001405.39_GRCh38.p13_genomic.gtf --sjdbOverhang 75 
```

## Then you need to align the reads to the genome
```
/home/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /mnt4/transciptome_compare/komarova/STAR_genome_index/ --runThreadN 10 --readFilesIn /mnt4/transciptome_compare/komarova/ERR3671306.fastq --outFileNamePrefix /mnt4/transciptome_compare/komarova/STAR_alignment/3671306 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```
## Counts table
After STAR, counts are obtained in files with .tab resolution and aligned reads in .out.bam. But it is better to run featureCount to calculate counts.
```
/home/tools/subread-2.0.1-source/bin/featureCounts -T 4 -s 0 -a /mnt4/transciptome_compare/komarova/GCF_000001405.39_GRCh38.p13_genomic.gtf -o /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts.txt /mnt4/transciptome_compare/komarova/STAR_alignment/*.out.bam
```
Where is:
- -T number of cores
- -s these data are "reverse"ly stranded
- -a path to annotation
- -o output file. The output of this tool is 2 files, a count matrix and a summary file that tabulates how many the reads were “assigned” or counted and the reason they remained “unassigned”.

## "Cleaning" of output file
Unnecessary columns need to be removed.
```
cut -f1,7-18 /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts.txt > /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts_Rmatrix.txt
```

## DESeq2 analysis
Differentially expressed gene analysis and signal pathway analysis was performed in RStudio. See HF_IB.Rmd script.
Differentially expressed genes were determined using a DESeq package in R with pairwise comparison of conditions. Only genes with Log2FoldChange>1 and padj<0.05were considered as differentially expressed. To find significantly enriched pathways we performed GO ontology analysis(gseaGO function), fast gene set enrichment analysis (fgsea function) and enrichment analysis (enrichGO function) using RStudio. 


## Results
We reprocessed RNA-seq data from the article and got just two differentially expressed genes (DEG). 

![Figure 1](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/Volcanoplot_HD_HF.png)

Due to the small number of DEG, only analysis of the functional groups of genes was performed. We used the gseGO and fgsea functions from clusterProfiler library in R. 
![Figure_2](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/GO_1.png)
In both cases, we got the same up and down regulated signal pathways in patients with heart failure. 
![Figure_3](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/FGSEA_1.png)
![Figure_4](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/FGSEA_2.png)
The pathways that control cellular respiration and protein localisation were downregulated in patients. 
Comparison of healthy donors and patients from our laboratory with fgsea gave many DEGs and slightly different corresponding signaling pathways. This can be explained by the not very good choice of the control group in the article and by another muscle.

PCA graph shows that the patients in the article are divided into two groups according to the first component. 
![Figure_5](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/PCA_plot_HD_HF.png)
![Figure_6](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/PCA_plot.png)
We have found that this is influenced by the left ventricular ejection fraction. The second group (HF_2) is characterized by a lower (almost critical) fraction compared to the first group (HF_1). In this case we also did GO analysis, FGSEA analysis and enrichment analysis on DEG. 
We identified more differentially expressed genes when comparing healthy donors and patients. 
![Figure_7](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/Volcanoplot_HF_HF.png)
Using GO analysis and FGSEA we found that left ventricular ejection fraction has an effect on RNA splicing, various metabolic pathways, immune response, cell adhesion and muscle contraction. Most of the signal pathways in HF_2 were downregulated. They included RNA splicing and cilium organization. 
![Figure_8](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/GO_2.png)
![Figure_9](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/FGSEA1_HF.png)
![Figure_10](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/FGSEA_2_HF.png)
Also, we have performed enrichment analysis on down and up regulated DEG. 
For DEG with Log2FoldChange<0:
![Figure_11](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/EnrichGO_HF_down.png)
For DEG with Log2FoldChange<0:
![Figure_12](https://github.com/Rita1612-GitHub/Transcriptomic_analysis_of_patients_with_heart_failure/blob/main/Pictures/EnrichGO_HF_up.png)

Thus, in the present work we have shown that the left ventricular ejection fraction has a key influence on the development of skeletal muscle atrophy in heart failure. 
