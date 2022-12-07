# analysis_BCR_DR
  
## Scripts for B cell receptor repertoire analysis
  
Immunosenescence is the decline in adaptive immune function with age that increases vulnerability to infection and participates in the feedback cycle that decreases health in old age.  
  
In this study we have 40 mice corresponding to 8 biological conditions. We performed BCR-Seq based on the work by Turchainova et al 2017 and prepared libraries from frozen spleen and ileum.  
  
```
General processing steps in analysis_BCR_DR/scripts/CMD_SLURM/ call into CMD_FOF:  
  
Unzip fastq files  
Filter by quality 20  
Map PCR2 primers (internal barcodes) to identify which reads correspond to each mouse  
Pair Forward and Reverse reads with the internal barcode information  
Split reads into fastqs corresponding to each mouse  
Extract UMI sequence from begining of Forward read  
Write the UMI information on the Reverse read  
Parse read headers and include information of file of origin  
Concatenate Forward and Reverse reads of the same mouse that come from different origins  
Build consensus from reads sharing the same UMI  
Pair consensus from forward and reverse  
Assemble forward and reverse into a merged full read  
Map PCR1 primers to identify the isotype of each molecule  
Transform to fasta  
Parse the isotype table  
Assign V-D-J genes  
Generate ChangeO database for analysis  
Extract only productive V-D-J rearrangements  
Use machine learning and hierarchical clustering to identify the threshold that splits each sample into clone and non-clone molecules  
Identify novel alleles and reassign genotypes including this information  
Using genotypes and threshold, define which reads correspond to each clone  
Create final germlines for analysis  
```
  
## Scripts for BCR data analysis  
We calculate and analyze many different parameters, gross and individually per isotype (and organ), and then all together we evaluate them in the context of the pathologies observed in the mice in study.  
  
\[SPL\|ILE\]\_parser\_tigger\_all.R = Find novel alleles  
\[SPL\|ILE\]\_novelAlleles\_plot.R = Plotting alleles identified during germline anotation for validation  
\[SPL\|ILE\]\_parser\_novelAlleles\_all.R = Plot all novel alleles and genotypes  
explore\_general.ipynb = Quality control, compare clonal assignment methods and thresholds, number of new alleles per sample  
hill\_diversity\_\[SPL\|ILE\]\_studio.Rmd = Estimate alpha diversity and clonal abundance  
\[SPL|ILE\]\_Alpha\_diversity.ipynb = Plotting Hill, Alpha boxplots, Kruskal with Mann Whitney, 2-way ANOVA  
isotypes\_alpha\_calculate.R = Estimate alpha diversity per isotype  
Isotype\_alpha\_\[SPL\|ILE\].ipynb = Plotting Hill, Alpha boxplots, 2-way ANOVA and linear regression through age, Kruskal with Mann Whitney, all by isotypes  
rdi\_beta\_\[SPL\|ILE\].Rmd = Use RDI to calculate beta diversity, whole and by isotype. PCOA and dendogram  
RDI\_plot\_analysis.ipynb = Plot boxplots, 2-way ANOVA and linear regression through age, Kruskall with Mann Whitney  
RDI\_isotypes.ipynb = Plot boxplots, 2-way ANOVA and linear regression through age, Anova between diets. By isotype  
Clonal\_Abundance.ipynb  = Plot log abundance curves, extract p20 (sum of frequencies from all clones with rank below or equal to 20), plot boxplots, Kruskall with Mann Whitney, 2-way ANOVA (no age effect, therefore no linear regression)  
ClonalAbundance\_isotypes\_calc.R = Calculate clonal abundance per isotype  
ClonalAbundance\_isotypes.ipynb = Plot log abundance curves, extract p20, boxplots, Kruskall with Mann Whitney, 2-way ANOVA (no age effect, therefore no linear regression)  
stacked\_isotypes.ipynb = Calculate frequency of isotypes, boxplots, stacked plots, Kruskall with Mann Whitney, 2-way ANOVA (no differences by age, so no linear regression), all by isotype  
mutational\_load\_quantification.R = Calculate mutational load for IMGT\_V for both synonymous and non-synonymous  
SHM\_study.ipynb = Stats for SHM calculated in R, 2-way ANOVA, Kruskall with Mann Whitney, boxplots, linear regression. Also for isotypes  
NaiveClassSwitch.ipynb = Find Class switch using Isotype and SHM status, linear regression, 2-way ANOVA, Kruskall with Mann Whitney  
CDR3\_length.ipynb = Plot gaussian curves by biorep, compare mu and sigma between groups using Kruskall with Mann Whitney, use KW to compare the full curves, compare ratio of number of significantly different biological replicates comparisons between ages and diets using linear regression and fisher tests  
CDR3\_length-isotypes.ipynb = Plot gaussian curves by biorep, compare mu and sigma between groups using Kruskall with Mann whithey, KW to compare full curves, ratio using linear regression and fishers, for all isotypes  
VJ\_geneUsage.Rmd = Calculate gene usage for tissues and isotypes  
VJusage.ipynb = Plot heatmaps for biorep and means, 2-way ANOVA  
VJusage-isotypes.ipynb = Plot heatmaps for biorep and means, Kruskall with Mann Whitney, linear regression, 2-way ANOVA for all isotypes  
Selection\_pressure.R = **testing**  
merging\_param\_patholo.ipynb = Formats F2 tumor load, SPL alpha diversity (richness, shannon, simpson), clonal abundance, SHM frequency (synonymous, non-synonymous, and mutational status 1 mutated, more than 0.01 mutation frequency, and 0 non mutated), isotypes relative frequency, class switch status, F2 pathology **missing many other things to calculate like RDI, isotype specific data, VJ usage etc**  
pathology\_alpha\_stats.Rmd = GLME of alpha diversity to pathology burden and presence/absence  
tumour\_alpha\_stats.Rmd = GLME of alpha diversity to tumor burden and presence/absence  
AlphaSPL\_pathology.ipynb = Boxplots of alpha diversity with specific pathologies/tumors to visualize association  
mixOmics\_igseqmicrobiome.Rmd = Diablo package to merge microbiome and spleen igseq clones, splsda, PC, auroc, circos, variables contributions, PC contributions  
mixomics\_ilesplmicro.Rmd = Mixomics also including ileum igseq, comparing clones with spleen and microbiome  
morbidity\_BCRparam.ipynb = Calculating morbidity score, OLS to p20, mufreq, mutstatus, relisotypes, naiveperc, CSperc, alpha; using rsquared for variance explained, spearman to know positive and negative correlation **This plots and calculations should be repeated with all the parameters to predict morbidity score**  

