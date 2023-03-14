---
title: "Analyses of APA dynamics accross human line"
---------------------------------------------------
author: "Mostafa Abouzaid"
--------------------------
date: "Last modified 2020-07-25"
--------------------------------

---

```{r warning=FALSE, include=FALSE}
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  echo = FALSE,
  comment = "#>"
)
```

# Overview

We used movAPA R package on a PAC-seq datasets generated from 8 Homosapien cell line from 3'end sequencing.
Here the poly(A) sites are already poly(A) site clusters (PACs) which were grouped from nearby cleavage sites.

# Preparations

## The PAC data

movAPA implemented the *PACdataset* object for storing the expression levels and annotation of PACs from various conditions/samples.
Almost all analyses of poly(A) site data in movAPA are based on the *PACdataset*.
The "counts" matrix is the first element in the array list of *PACdataset*, which stores non-negative values representing expression levels of PACs.
The "colData" matrix records the sample information and the "anno" matrix stores the genome annotation or additional information of the poly(A) site data.
for removing internal priming or poly(A) signal analyses.
we used reference genome sequences that are represented as a *BSgenome* object or stored in a fasta file.The function *parseGff* is used to parse the given annotation.
The function *annotatePAC* is used to annotate a *PACdataset* with a GFF/GTF file.

```{r include=FALSE, eval=TRUE}
setwd("J:/SchafferLab/wagner/data/cleandata/human/decomp/count")
library(movAPA)
library(kableExtra)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(ComplexHeatmap, quietly = TRUE)
library(tibble)
library('biomaRt')

theme=theme(axis.text.x = element_text(angle = 90, hjust = 1))
ht_opt$message = FALSE
library("BSgenome.Hsapiens.UCSC.hg38")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
TxDb <-TxDb.Hsapiens.UCSC.hg38.knownGene
gffFile="Homo_sapien.chr.gtf"
gff=parseGff(gffFile) 
source("C:\\Users\\Mosta\\GitHub\\movAPA\\R\\R_funclib_movAPA.r")
#gff=parseTxdb(TxDb)
```

```{r, eval=TRUE}
cond1=c('wt1','L352','het3','het4')
cond2=c('mut1',"L412",'mut3','mut4')
pac=read.table("all.PAC.PAT_chr_count", header=F, stringsAsFactors = F)
header=read.table("all.PAC.header", header=T, stringsAsFactors = F)
colnames(pac)=colnames(header)
colData=read.delim("col.data",row.names = 1)
row_sub = apply(pac[,c(9:24)], 1, function(row) all(row <= 100 , row >= 1))
pac2=pac[row_sub,]
##Subset as usual
#colnames(pac)=c('chr','coord','x','strand')
#pac=pac[,c('chr','strand','coord')]
PACds=readPACds(pacFile=pac2, colDataFile=colData, noIntergenic=TRUE,PAname = 'PA' )


```

# Preprocessing of PAC data

Internal priming (IP) artifacts can be removed by the *removePACdsIP* function.
Here, PACs with six consecutive or more than six As within the -10 to +10 nt window are considered as internal priming artifacts.
We scaned the internal priming artifacts in PACds and get two *PACdatasets* recording internal priming PACs and real PACs.
while The function *mergePACds* used to group nearby cleavage sites into PACs.
*annotatePAC is used* to annotate a *PACdataset* with a GFF/GTF file .

```{r fig.dim = c(6, 4), fig.align='left'}
PACds=annotatePAC(PACds,gff)
PACdsIP=removePACdsIP(PACds, bsgenome, returnBoth=TRUE, 
                      up=-10, dn=10, conA=6, sepA=7)

IP=data.frame(per=c(length(PACdsIP$real)/(length(PACdsIP$real)+length(PACdsIP$ip))*100,length(PACdsIP$ip)/(length(PACdsIP$real)+length(PACdsIP$ip))*100),PACS=c("Real Pacs (RP)","Internal priming (IP)"))
 
```

```{r results='hide'}

ggpubr::ggbarplot(data = IP,y = "per",x="PACS",fill ="PACS",xlab = "",ylab = "(%" ,sort.val = "desc",rotate=TRUE)

```

## Normalization

The function *normalizePACds* called for normalization, which implements three strategies including TPM (Tags Per Million), the normalization method of DESeq ([Anders and Huber, 2010](#3)), and the TMM method used in EdgeR ([Robinson, et al., 2010](#4)).

```{r fig.align='left', fig.height=6, fig.width=4, message=FALSE, warning=FALSE, results='hide'}
## Here normalization method TMM (or EdgeR) is used, 
## while you may also choose TPM or DESeq.
#PACdsClust=PACds
par(family="serif",face="bold", size=1)
librarySizes <- data.frame(librarySizes=colSums(PACds@counts),
                           sample=names(colSums(PACds@counts)))
lib=cbind(colData,librarySizes)
PACds=normalizePACds(PACds, method='TPM')
counts=as.matrix(PACds@counts)
logcounts <- as.data.frame(vst(counts,nsub = 50))
col=colData %>% rownames_to_column("variable") 
box=logcounts %>% melt(variable_name = "Sample") %>% left_join(col)

ggbarplot(lib,x ="sample",
          y="librarySizes",
          fill = "group",
          short.panel.labs = TRUE,
          sort.val = "desc") + theme(axis.text.x = element_text(angle = 90)) +facet_wrap(~group,scales = "free")


ggboxplot(box,x="variable",y="value",
        xlab="",
        bxp.errorbar = TRUE,
        notch = TRUE,
        fill = "group",
        short.panel.labs = TRUE,
        ylab="Log2(Counts)")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+facet_wrap(~group,scales = "free") 
```

# Principle Component Analysis

A principle components analysis (PCA) is an example of an unsupervised analysis, where we don't specify the grouping of the samples.
If your experiment is well controlled and has worked well, we should that replicate samples cluster closely, whilst the greatest sources of variation in the data should be between treatments/sample groups.
It is also an incredibly useful tool for checking for outliers and batch effects.

To run the PCA we should first normalise our data for library size and transform to a log scale.DESeq2 provides two commands that can be used to do this, here we will use the command `rlog`.
`rlog` performs a log2 scale transformation in a way that compensates for differences between samples for genes with low read count and also normalizes between samples for library size.

You can read more about `rlog`, it's alternative `vst` and the comparison between the two [here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations).

To plot the PCA results we will use the `autoplot` function from the `ggfortify` package [@Tang2016].
`ggfortify` is built on top of `ggplot2` and is able to recognise common statistical objects such as PCA results or linear model results and automatically generate summary plot of the results in an appropriate manner.

```{r pcaPlot, fig.align="center", message=FALSE, warning=FALSE}


# run PCA
PCA <-  prcomp(t(logcounts))

autoplot(PCA,
         data = col, 
         colour="group", 
         shape=10,
         label = TRUE,
         labe_size = 1,
         frame = TRUE, 
         frame.type = 'norm',
         size=1)+theme_minimal()
minkowski=hclust(dist(t(logcounts),method = "minkowski"))
plot(minkowski,main="distance",xlab="minkowski")
```

# Statistical analyses of PACs

To make statistics of distributions of PACs for each sample, first we pooled replicates.Then we calculated statistics of distribution of PACs using different PAT cutoffs.
minPAT=5 means that only PACs with \>=5 reads are used for statistics.

```{r}
PACds1=subsetPACds(PACds, group='group', pool=TRUE)
pstats=movStat(PACds1, minPAT=c(1, 5, 10, 20, 50, 60), ofilePrefix=NULL)

pstats$pat10 %>%
  kbl(caption = "Table 1") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Statistical results can be visualized by barplots to show PAC#, PAT#, APA gene%, PAC%, PAT% across samples and genomic regions.
Here we plot all statistical results with cutoffs 5 and 10, with each plot having two smaller plots corresponding to the two cutoffs.

```{r }

plotPACdsStat(pstats, pdfFile=NULL, minPAT=c(5,10)) 
```

Plot specific cutoffs and conditions.

```{r eval=FALSE}
plotPACdsStat(pstats, pdfFile=NULL, 
              minPAT=c(5,10), conds=c(cond1,cond2))
```

Plot the overall distributions using pooled samples (total) and two cutoffs.

```{r }
plotPACdsStat(pstats, pdfFile=NULL, 
              minPAT=c(5,10), conds=c('total'))
```

Plot the overall distributions using pooled samples (total) and one cutoff.

```{r }
plotPACdsStat(pstats, pdfFile=NULL, 
              minPAT=c(10), conds=c('total'))
```

# Poly(A) signals and sequences

movAPA provides several functions, including *annotateByPAS*, *faFromPACds*, *kcount*, and *plotATCGforFAfile*, for sequence extraction and poly(A) signal identification.

## Poly(A) signals

Annotate PACs by corresponding signal of AATAAA located upstream 50 bp of the PAC.

```{r}
PACdsPAS=annotateByPAS(PACds, bsgenome, grams='AATAAA', 
                       from=-50, to=-1, label=NULL)

```

Scan AATAAA's 1nt variants.

```{r}
PACdsPAS=annotateByPAS(PACds, bsgenome, grams='V1', 
                       from=-50, to=-1, label=NULL)
table(PACdsPAS@anno$V1_gram)  %>%
  kbl(caption = "Table 2") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Scan custom k-grams.

```{r}
PACdsPAS=annotateByPAS(PACds, bsgenome, 
                       grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
                       from=-50, to=-1, label='GRAM')
table(PACdsPAS@anno$GRAM_gram)  %>%
  kbl(caption = "Table 3") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Scan patterns with priority: AATAAA \>\> ATTAAA \>\> remaining k-grams.

```{r}
PACdsPAS=annotateByPAS(PACds, bsgenome,
                       grams=c('AATAAA','ATTAAA','GATAAA','AAAA'),
                       priority=c(1,2,3,3), 
                       from=-50, to=-1, label='GRAM')
table(PACdsPAS@anno$GRAM_gram)  %>%
  kbl(caption = "Table 4") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Plot signal logos.

```{r message=FALSE, warning=FALSE}
pas=PACdsPAS@anno$GRAM_gram[!is.na(PACdsPAS@anno$GRAM_gram)]
plotSeqLogo(pas)
```

Here we show another example to scan mouse signals in human PACs.
First, we get mouse signals and set the priority.

```{r}
v=getVarGrams('mm')
priority=c(1,2,rep(3, length(v)-2))
```

Then scan upstream regions of PACs for mouse signals.

```{r}
PACdsMM=annotateByPAS(PACds, bsgenome, grams=v, 
                      priority=priority, 
                      from=-50, to=-1, label='mm')
```

Prepare the data to plot PAS distributions.

```{r results='hide'}
pas=as.data.frame(cbind(region=PACdsMM@anno$ftr, PAS=PACdsMM@anno$mm_gram))
pas$PAS[is.na(pas$PAS)]='NOPAS'
pas$PAS[pas$PAS %in% v[-c(1:2)]]='Variants'
n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
n=merge(n, n2)
n$PAC=n$nPAC/n$nTot
n=n[n$PAS!='NOPAS', ]
n$PAS=factor(n$PAS, levels=rev(c('AATAAA', 'ATTAAA','Variants', 'NOPAS')))
n$region=factor(n$region, 
                levels=c('3UTR','Ext_3UTR', 'intergenic','intron','CDS','5UTR'))
```

Plot PAS distributions.

```{r fig.dim = c(6, 4), fig.align='left'}
ggplot(data=n, aes(x=region, y=PAC, fill=PAS)) + 
  geom_bar(stat="identity") + 
  ylab("PAC Fraction") + theme_bw()
```

## Exract sequences

The *faFromPACds* function used various options to extract sequences of interest.

```{r eval=FALSE}

## Extract the sequence of PACs, from UPA_start to UPA_end. 
faFromPACds(PACds, bsgenome, what='pac', fapre='pac')

## Extract upstream 300 bp ~ downstream 100 bp around PACs, 
## where the position of PAC is 301.
faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
            up=-300, dn=100)

## Divide PACs into groups of genomic regions and then extract sequences for each group.
faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
            up=-100, dn=100, byGrp='ftr')

## Extract sequences for only 3UTR PACs.
faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
            up=-300, dn=100, byGrp=list(ftr='3UTR'))

## Extract sequences for only 3UTR PACs and separate sequences by strand.
faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
            up=-300, dn=100,
            byGrp=list(ftr='3UTR', strand=c('+','-')))

faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
                    up=-300, dn=100, byGrp='ftr')

## Extract sequences of genomic regions where PACs are located.
faFromPACds(PACds, bsgenome, what='region', fapre='region', byGrp='ftr')
faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn', 
                    up=-300, dn=100, byGrp='ftr')
```

## Base compostions and k-grams

The function *plotATCGforFAfile* is for plotting single nucleotide profiles for given fasta file(s), which is particularly useful for discerning base compositions surrounding PACs.
First trim sequences surrounding PACs.
Sequences surrounding PACs in different genomic regions are extracted into files.
The PAC position is 301.Then plot base compositions for specific sequence file(s).

```{r fig.dim=c(6,4)}
faFiles=c("updn.3UTR.fa", "updn.5UTR.fa","updn.exon.fa","updn.CDS.fa", "updn.intergenic.fa", "updn.intron.fa")
## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
plotATCGforFAfile (faFiles, ofreq=FALSE, opdf=FALSE, 
                   refPos=301, mergePlots = TRUE)
```

After extracting sequences, we can call the *kcount* function to obtain the number of occurrences or frequencies of k-grams from the whole sequences or a specified region of sequences.
Particularly, specific k-grams (e.g., AAUAAA, AUUAAA) or a value of k (e.g., k=6 means all hexamers) can be set.

```{r}
## Count top 10 hexamers (k=6) in the NUE region 
## (normally from 265~295 if the PAC is at 301).
fafile='updn.3UTR.fa'
kcount(fafile=fafile, k=6, from=265, to=295, topn=10)  %>%
  kbl(caption = "Table 5") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## Count given hexamers.
kcount(fafile=fafile, grams=c('AATAAA','ATTAAA'), 
       from=265, to=295, sort=FALSE)  %>%
  kbl(caption = "Table 6") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## Count AATAAA and its 1nt variants in a given region.
kcount(fafile=fafile, grams='v1', from=265, to=295, sort=FALSE)  %>%
  kbl(caption = "Table 7") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

# Quantification of PACs by various metrics

movAPA provides various metrics to measure the usages of PACs across samples, including three metrics for the quantification of the usage of each single poly(A) site by the *movPAindex* function and four metrics for the quantification of APA site usage of a gene by the *movAPAindex* function.

## Quantification of each PAC by *movPAindex*

*movPAindex* provides three metrics for the quantification of each PAC in a gene, including "ratio", "Shannon", and "geo".
First you can merge replicates of the same sample and remove lowly expressed PACs before calculate the index.

```{r message=FALSE}
p=subsetPACds(PACds, group='group', pool=TRUE, totPACtag=20)
```

the line-specificity.
Q or H=0 means that the PAC is only expressed in one line.
NA means the PAC is not expressed in the respective line.Use the relative expression levels (ratio) to calculate line-specificity.

```{r}
paShan=movPAindex(p, method='shan', shan.ratio = TRUE) 
head(paShan, n=10)  %>%
  kbl(caption = "Table 10") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Cacluate the geo metric, which is only suitable for APA genes.
NA means no PAC of the gene is expressed in the respective tissue.
geo\>0 means the PAC is used more than average usage of all PACs in the gene.
geo\~0 means similar usage; \<0 means less usage.

```{r}
paGeo=movPAindex(p, method='geo')
head(paGeo, n=10)  %>%
  kbl(caption = "Table 11") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Cacluate the ratio metric, which is only suitable for APA genes.
NA means no PAC of the gene is expressed in the respective tissue.

```{r}
paRatio=movPAindex(p, method='ratio')
head(paRatio)  %>%
  kbl(caption = "Table 12") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Plot a heatmap to show the distribution of line-specificity of PACs.
It is only reasonable to plot the heatmap of the Shanno metric.
Or you may filter the proximal or distal PAC of the gene first and plot the ratio or geo metrics.
First, remove rows with NA and then plot the heatmap.

```{r fig.dim=c(6,6), message=FALSE, warning=FALSE}
paShanHm=paShan[, -(1:3)]
paShanHm=scale(paShanHm[rowSums(is.na(paShanHm))==0, ])
Heatmap(paShanHm, show_row_names=FALSE,  cluster_columns = FALSE, 
        heatmap_legend_param = list(title = 'line\nspecificity'))
```

Calculate the Genotype-specificity for each replicate.

```{r}
paShan=movPAindex(PACds, method='shan',shan.ratio = TRUE)
```

```{r fig.dim=c(6,6)}
## Plot heamap to show the consistency among replicates.
paShanHm=paShan[, -(1:3)]
paShanHm=scale(paShanHm[rowSums(is.na(paShanHm))==0, ])
Heatmap(paShanHm, show_row_names=FALSE,  cluster_columns = TRUE, 
        heatmap_legend_param = list(title = 'line\nspecificity'))  
```

## Quantification of APA by *movAPAindex*

The *movAPAindex* function provides four gene-level metrics for the quantification of APA site usage, including RUD (Relative Usage of Distal PAC) ([Ji, et al., 2009](#5)), WUL (Weighted 3' UTR Length) ([Ulitsky, et al., 2012](#6); [Fu, et al., 2016](#1)), SLR (Short to Long Ratio) ([Begik, et al., 2017](#7)), and GPI (Geometric Proximal Index) ([Shulman and Elkon, 2019](#8)).

Get APA index using the RUD method.

```{r eval=FALSE}
geneRUD=movAPAindex(p, method="RUD",
                    choose2PA=NULL, RUD.includeNon3UTR=TRUE)
geneRUD=scale(geneRUD[rowSums(is.na(geneRUD))==0, ])
head(geneRUD, n=10)%>%
  kbl(caption = "Table 13") %>%
  kable_classic(full_width = F, html_font = "Cambria")
Heatmap(geneRUD, show_row_names=FALSE,  cluster_columns = F, 
        heatmap_legend_param = list(title = 'RUD'))
```

Get APA index by method=SLR, using the proximal and distal PACs.

```{r eval=FALSE}
geneSLR=movAPAindex(p, method="SLR", choose2PA='PD',RUD.includeNon3UTR = FALSE)
geneSLR=scale(geneSLR[rowSums(is.na(geneSLR))==0, ])
head(geneSLR, n=10) %>%
  kbl(caption = "Table 14") %>%
  kable_classic(full_width = F, html_font = "Cambria")
suppressMessages(
Heatmap(geneSLR, show_row_names=FALSE))
```

Get APA index by method=GPI, using the proximal and distal PACs.

```{r fig.dim=c(6,6), warning=FALSE}
geneGPI=movAPAindex(PACds, method="GPI", choose2PA='PD')
head(geneGPI) %>%
  kbl(caption = "Table 15") %>%
  kable_classic(full_width = F, html_font = "Cambria")
geneGPI=geneGPI[rowSums(is.na(geneGPI))==0, ]
Heatmap(geneGPI, show_row_names=FALSE,  cluster_columns = TRUE, 
        heatmap_legend_param = list(title = 'GPI'))
```

# DE genes

3' seq data have been demonstrated informative in quantifying expression levels of genes by summing up 3' seq reads of all PACs in a gene ([Lianoglou, et al., 2013](#9)).
To detect DE genes between samples with 3' seq, we implemented the function *movDEgene* with the widely used R package DESeq2.

## Detect DE genes

First we show an example of detecting DE genes for two conditions.

```{r message=FALSE, warning=FALSE, results='hide'}
## Subset two conditions first.
PACds=subsetPACds(PACds, group='group',cond1=cond1, cond2=cond2)
## Detect DE genes using DESeq2 method, 
## only genes with total read counts in all samples >=50 are used.
DEgene=movDEGene(PACds=PACds, method='DESeq2', group='group', minSumPAT=10)
```

Make statistics of the DE gene results; genes with padj\<0.05 & log2FC\>=0.5 are considered as DE genes.

```{r}
stat=movStat(object=DEgene, padjThd=0.05, valueThd=0.5)
stat$nsig %>%
  kbl(caption = "Table 16") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

We can also detect DE genes among more than two conditions.

```{r}
## Number of DE genes in each pair of conditions.
stat$nsig %>%
  kbl(caption = "Table 17") %>%
  kable_classic(full_width = F, html_font = "Cambria")
## Overlap between condition pairs.
stat$ovp %>%
  kbl(caption = "Table 18") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

## Output DE genes

Output *movStat* results into files: "DEgene.plots.pdf" and 'DEgene.stat'.
Several heatmaps are generated.

```{r eval=FALSE}
outputHeatStat(heatStats=stat, ostatfile='DEgene.stat', plotPre='DEgene')
```

You can further call movSelect() to select DE gene results with more information.
Select DE gene results with full information including the read counts in each sample.

```{r message=FALSE}
selFull=movSelect(DEgene, condpair='wt1.mut1', padjThd=0.05, valueThd=1, 
                  out='full', PACds=PACds)

```

Select DE gene results with only padj and value.
Here value is log2.

```{r}
sel=movSelect(DEgene, condpair='wt1.mut1', 
              padjThd=0.05, valueThd=1, out='pv')
head(sel)  %>%
  kbl(caption = "Table 19") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Output gene names of DE genes.

```{r}
sel=movSelect(DEgene, condpair='wt1.mut1', 
              padjThd=0.05, upThd=0.5, out='gene')

```

# DE PACs

movAPA provides the function *movDEPAC* to identify DE PACs between samples.
Three strategies were utilized: (i) using DESeq2 with replicates; (ii) using DEXseq with replicates; (iii) using chi-squared test without replicates ("chisq").

## Detect DE PACs

First we show an example of detecting DE PACs among all pairwise conditions using three different methods.
Only PACs with total read counts in all samples \>=20 are used.

```{r results='hide', message=FALSE}
DEPAC=movDEPAC(PACds, method='DESeq2', group='group', minSumPAT=20)
DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=20)
DEqPAC=movDEPAC(PACds, method='chisq', group='group', minSumPAT=20)
```

Number of DE PACs among methods.

```{r results='hide'}
## Get significant DE results.
stat1=movStat(object=DEPAC, padjThd=0.05, valueThd=1)
stat2=movStat(object=DEXPAC, padjThd=0.05, valueThd=1)
stat3=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)
```

```{r fig.dim=, warning=FALSE}
## Count the number of DE PACs by different methods.
nsig=as.data.frame(cbind(stat1$nsig, stat2$nsig, stat3$nsig))
colnames(nsig)=c('DESeq2','DEXseq','Chisq.test')
nsig$contrast=rownames(nsig)

## Plot a barplot.
nsig=melt(nsig, variable.name='Method')
nsig %>%
  kbl(caption = "Table 20") %>%
  kable_classic(full_width = F, html_font = "Cambria")
ggplot(data=nsig, aes(x=contrast, y=value, fill=Method)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("DE PAC#") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

We can also detect DE PACs between two given conditions.

```{r eval=FALSE}
## First subset PACs in two conditions.
PACds1=subsetPACds(PACds, group='group', 
                   cond1='mut1', cond2='wt1', choosePA='apa')
## Detect DE PACs.
DEPAC1=movDEPAC(PACds1, method='DESeq2', group='group', minSumPAT=10)
DEXPAC1=movDEPAC(PACds1, method='DEXseq', group='group', minSumPAT=10)
DEqPAC1=movDEPAC(PACds1, method='chisq', group='group', minSumPAT=10)
```

## Statistics of DE PACs

Make statistics of the DE PACs result by DESeq2 method (*DEPAC*).

```{r results='hide'}
stat=movStat(object=DEPAC, padjThd=0.05, valueThd=1)
```

```{r}
## Number of DE PACs between conditions.
stat$nsig  %>%
  kbl(caption = "Table 20") %>%
  kable_classic(full_width = F, html_font = "Cambria")
## Overlap of DE PACs between different pairs of conditions.
head(stat$ovp) %>%
  kbl(caption = "Table 21") %>%
  kable_classic(full_width = F, html_font = "Cambria")

```

Stat the DE PAC result from the chisq-test method, here the value column of DEqPAC is 1-pvalue_of_the_gene.
So using padjThd=0.05 and valueThd=0.95 means filtering DE PACs with adjusted pvalue of PAC \<0.05 and adjusted pvalue of gene \<0.05.

```{r eval=FALSE}
stat=movStat(object=DEqPAC, padjThd=0.05, valueThd=0.95)
```

## Output DE PACs

We can use *movSelect* to output full or simple list of DE PACs.

```{r message=FALSE, warning=FALSE}
## Here method is DEXseq, so the valueThd (log2FC) threshold is automatelly determined. 
sel=movSelect(aMovRes=DEXPAC, condpair='wt1.mut1', 
              padjThd=0.05, out='full', PACds=PACds)
head(sel, n=2)%>%
  kbl(caption = "Table 22") %>%
  kable_classic(full_width = F, html_font = "Cambria")

sel=movSelect(aMovRes=DEPAC, condpair='wt1.mut1', 
              padjThd=0.05, out='full', PACds=PACds)
head(sel, n=2)%>%
  kbl(caption = "Table 23") %>%
  kable_classic(full_width = F, html_font = "Cambria")


sel=movSelect(aMovRes=DEqPAC, condpair='wt1.mut1', 
              padjThd=0.05, out='full', PACds=PACds)
head(sel, n=2)%>%
  kbl(caption = "Table 24") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

## Visualize DE PACs in a gene

Here we take the DEPAC result for example to show the visualization of DE PACs in a gene.
Visualize DE PACs in an example gene by *movViz*.
First, we examine all PACs in this gene.
There are three PACs, two in 3'UTR and one in extended 3'UTR.
But the expression level of the PAC in extended 3UTR is only 3.

```{r}
gene='ENST00000304858'
gp=PACds[PACds@anno$gene==gene, ]
gp=
cbind(gp@anno$ftr, rowSums(gp@counts),make.row.names = FALSE)
```

Visualize PACs of this gene in individual conditions.
Here the Y-axis is read count, the scale of which is different among conditions.DE PACs identified by DESeq2 method with padj \< padjThd are highlighted in dashed yellow lines.

```{r fig.height=4, fig.width=6, message=FALSE, warning=FALSE}
outputHeatStat(heatStats=stat, ostatfile='DEPAC.stat', plotPre='DEPAC', 
               show_rownames = TRUE)
movViz(object=DEPAC, gene=gene, txdb=, PACds=PACds, collapseConds=FALSE, 
       padjThd=0.01, showRatio=FALSE, showAllPA=TRUE) 
```

We can also show condition pairs in individual tracks and only display and/or highlight given condition pairs.
If padjThd is given, then the DE PACs (padj \< padjThd) will be highlighted (dashed yellow line).

```{r fig.dim=c(7, 5), message=FALSE, warning=FALSE}
cond=paste(cond1,cond2,sep = ".")
movViz(object=DEPAC, gene=gene, txdb=gff, PACds=PACds, collapseConds=TRUE, 
       padjThd=0.01, showPV=TRUE, showAllPA=T, showRatio=T,
       conds=DEPAC@conds[cond, ], highlightConds=DEPAC@conds[cond, ])
```

# APA site switching

The function *movAPAswitch* is used to detect both canonical and non-canonical APA site switching events.
The strategy of *movAPAswitch* is similar to the strategy based on DE PACs in *movUTRtrend* but with higher flexibility.
If a gene has more than two PACs, then each pair of PACs (denoted as PA1 and PA2) are analyzed.
The following criteria are used to determine a APA switching event: whether PA1 or PA2 are DE; average read count for both sites; distance between PA1 and PA2; average read count for a gene; relative change of PA1 and PA2 (RC); read count ratio (PA1:PA2) \>1 in one sample and \<1 in another sample; p-value of the Fisher' s exact test for PA1 and PA2 read counts between samples.
Pairs of PACs that meet user specified conditions are considered as APA site switching events.
Users can use the *movSelect* function to filter 3' UTR switching events or APA site switching events with higher flexibility.

## Detect 3'UTR-PAC switching

First get DE PAC results by DEXseq.

```{r eval=FALSE}
DEXPAC=movDEPAC(PACds, method='DEXseq', group='group', minSumPAT=10)
```

Then get 3'UTR switching genes, usig selectOne=NULL to detect all pairs of switching PACs.

```{r message=FALSE, warning=FALSE, results='hide'}
swDEX=movAPAswitch(PACds, group='group',aMovDEPACRes=DEXPAC,
                   avgPACtag=5, avgGeneTag=10,
                   only3UTR=TRUE,
                   DEPAC.padjThd=0.1, nDEPAC=1,
                   mindist=50, fisherThd=0.1, logFCThd=0.5, 
                   cross=FALSE, selectOne=NULL)
```

Stat the switching results.

```{r message=FALSE}
stat=movStat(object=swDEX, padjThd=0.1, valueThd=1)
stat$nsig   %>%
  kbl(caption = "Table 25") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Output switching genes with full information

```{r message=FALSE}
sel=movSelect(aMovRes=swDEX, condpair='wt1.mut3', 
              padjThd=0.1, valueThd=1, out='full')
head(sel, n=2)  %>%
  kbl(caption = "Table 26") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

## Visualization of 3'UTR-PAC switching

Show one switching gene (ENST00000376802), where switching happens between a 3'UTR PAC a.
This gene has 2 PAC in 3UTR; the APA-site switching happens between mut1\~wt1

```{r}
gene='ENST00000376802'
gp=PACds[PACds@anno$gene==gene, ]
cbind(gp@anno$ftr, rowSums(gp@counts))%>%
  kbl(caption = "Table 27") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Plot all PACs of this gene in all conditions and replicates.
Highlight PACs involving in the switching analysis in orange.

```{r  fig.height=4, message=FALSE, warning=FALSE}
movViz(object=swDE, gene=gene, txdb=gff, PACds=PACds, 
       showRatio=TRUE, padjThd=0.01, showAllPA=FALSE)
```

## Detect APA-site switching

Detect APA switching events involving non-3'UTR PACs, using selectOne=NULL to get all pairs of switching PACs.

```{r message=FALSE, results='hide'}
swDE=movAPAswitch(PACds, group='group', aMovDEPACRes=DEXPAC,
                   avgPACtag=10, avgGeneTag=20,
                   only3UTR=FALSE,
                   DEPAC.padjThd=0.1, nDEPAC=1,
                   mindist=50, fisherThd=0.1, logFCThd=0.5, 
                  cross=FALSE, selectOne=NULL)
```

Stat the switching results.

```{r message=FALSE, warning=FALSE}
stat=movStat(object=swDE, padjThd=0.1, valueThd=1)
stat$nsig %>%
  kbl(caption = "Table 26") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Output switching genes with full information for anther\~embryo.

```{r message=FALSE}
sw=movSelect(aMovRes=swDE, condpair='wt1.mut1', 
             padjThd=0.01, valueThd=1, out='full')
head(sw[order(sw$fisherPV), ], n=10) %>%
  kbl(caption = "Table 27") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

## Subset PACds by switching genes or PACs

First get list of genes or PACs of switching events, then subset PACds by genes or PACs.

```{r message=FALSE}
genes=movSelect(aMovRes=swDE, condpair='wt1.mut1', 
                padjThd=0.01, valueThd=1, out='gene')
swPAC=subsetPACds(PACds, genes=genes, verbose=TRUE)
table(swPAC@anno$ftr)

PAs=movSelect(aMovRes=swDE, condpair='wt1.mut1', padjThd=0.01,
              valueThd=1, out='pa')
swPAC=subsetPACds(PACds, PAs=PAs, verbose=TRUE)
table(swPAC@anno$ftr) %>%
  kbl(caption = "Table 22") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

## Visualization of APA-site switching

Show one switching gene (ENST00000423485), where switching happens between a 3'UTR PAC and an exonic PAC.
This gene has 2 PACs in CDS and 1 PAC in 3UTR; the APA-site switching happens between wt1.mut1.

```{r}
gene='ENST00000423485'
gp=PACds[PACds@anno$gene==gene, ]
cbind(gp@anno$ftr, rowSums(gp@counts))%>%
  kbl(caption = "Table 23") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

Plot all PACs of this gene in all conditions and replicates.
Highlight PACs involving in the switching analysis in orange.

```{r  fig.height=4, message=FALSE, warning=FALSE}
movViz(object=swDE, gene=gene, txdb=gff, PACds=PACds, 
       showRatio=TRUE, padjThd=0.01, showAllPA=FALSE)
```

# Session Information

The session information records the versions of all the packages used in the generation of the present document.

```{r}
sessionInfo()
```

# References {#refer}

[[1] Fu, H., Yang, D., Su, W., et al. Genome-wide dynamics of alternative polyadenylation in rice. Genome Res. 2016;26(12):1753-1760.]{#1}

[[2] Zhu, S., Ye, W., Ye, L., et al. PlantAPAdb: A Comprehensive Database for Alternative Polyadenylation Sites in Plants. Plant Physiol. 2020;182(1):228-242.]{#2}

[[3] Anders, S. and Huber, W. Differential expression analysis for sequence count data. Genome Biol. 2010;11(10):2010-2011.]{#3}

[[4] Robinson, M.D., McCarthy, D.J. and Smyth, G.K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 2010;26(1):139-140.]{#4}

[[5] Ji, Z., Lee, J.Y., Pan, Z., et al. Progressive lengthening of 3' untranslated regions of mRNAs by alternative polyadenylation during mouse embryonic development. Proc. Natl. Acad. Sci. USA 2009;106(17):7028-7033.]{#5}

[[6] Ulitsky, I., Shkumatava, A., Jan, C.H., et al. Extensive alternative polyadenylation during zebrafish development. Genome Res. 2012;22(10):2054-2066.]{#6}

[[7] Begik, O., Oyken, M., Cinkilli Alican, T., et al. Alternative Polyadenylation Patterns for Novel Gene Discovery and Classification in Cancer. Neoplasia 2017;19(7):574-582.]{#7}

[[8] Shulman, E.D. and Elkon, R. Cell-type-specific analysis of alternative polyadenylation using single-cell transcriptomics data. Nucleic Acids Res 2019;47(19):10027-10039.]{#8}

[[9] Lianoglou, S., Garg, V., Yang, J.L., et al. Ubiquitously transcribed genes use alternative polyadenylation to achieve tissue-specific expression. Genes Dev. 2013;27(21):2380-2396.]{#9}
