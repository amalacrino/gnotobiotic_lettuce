# Data analysis

## Load libraries

```R
library("phyloseq")
library("ggplot2")
library("vegan")
library("ape")
library("dplyr")
library("picante")
library("emmeans")
library("car")
library("lme4")
library("data.table")
library("microbiome")
library("RVAideMemoire")
library("BiocParallel")
library("ggpubr")
library("phytools")
library("RColorBrewer")
library("decontam")
library("Maaslin2")
library("Wrench")
```

## Load data

```R
load(file = 'data/16s/ASV_table.rds')
ps.16s <- ps
track.16s <- track
load(file = 'data/its/ASV_table.rds')
ps.its <- ps
track.its <- track
metadata <- read.table("metadata.txt", sep = "\t", header = T, row.names = 1)
sample_data(ps.16s) <- metadata
sample_data(ps.its) <- metadata
tax_table(ps.its) <- gsub(".*_","",tax_table(ps.its))
```

## Remove contaminants and cleanup the dataset

```R
ps.16s <- subset_taxa(ps.16s, Class !="Chloroplast")
ps.16s <- subset_taxa(ps.16s, Order !="Chloroplast")
ps.16s <- subset_taxa(ps.16s, Family !="Mitochondria")

ps.its <- subset_taxa(ps.its, Order !="Chloroplast")
ps.its <- subset_taxa(ps.its, Family !="Mitochondria")

remove.cont <- function(ps){
  sample_data(ps)$is.neg <- sample_data(ps)$compartment == "inoc_NTC"
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.05)
  cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
  cont.remove <- row.names(cont.remove)
  allTaxa = taxa_names(ps)
  allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
  ps <-  prune_taxa(allTaxa, ps)
  temp <- sample_data(ps)[["group"]] %ni% "inoc_NTC"
  ps <- prune_samples(samples = temp, ps)
  return(ps)
}

ps.16s <- remove.cont(ps.16s)
ps.its <- remove.cont(ps.its)

temp <- sample_data(ps.16s)[["group"]] %in% c("inoc_46", "inoc_47", "inoc_48", "inoc_49")
ps.16s.ctrl <- prune_samples(samples = temp, ps.16s)
temp <- sample_data(ps.its)[["group"]] %in% c("inoc_46", "inoc_47", "inoc_48", "inoc_49")
ps.its.ctrl <- prune_samples(samples = temp, ps.its)

temp <- sample_data(ps.16s)[["group"]] %ni% c("inoc_46", "inoc_47", "inoc_48", "inoc_49")
ps.16s <- prune_samples(samples = temp, ps.16s)
temp <- sample_data(ps.its)[["group"]] %ni% c("inoc_46", "inoc_47", "inoc_48", "inoc_49")
ps.its <- prune_samples(samples = temp, ps.its)

ps.16s <- prune_samples(sample_sums(ps.16s)>=1000, ps.16s)
ps.its <- prune_samples(sample_sums(ps.its)>=1000, ps.its)
ps.16s <- filter_taxa(ps.16s, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.its <- filter_taxa(ps.its, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.16s.ctrl <- filter_taxa(ps.16s.ctrl, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps.its.ctrl <- filter_taxa(ps.its.ctrl, function (x) {sum(x > 0) > 1}, prune=TRUE)
```

## Normalize data

```R
wrench.norm.ps <- function(ps){
  count_tab <- as.matrix(t(data.frame(otu_table(ps))))
  group <- sample_data(ps)$compartment
  W <- wrench(count_tab, condition=group)
  norm_factors <- W$nf
  head(norm_factors)
  norm_counts <- sweep(count_tab, 2, norm_factors, FUN = '/')
  norm_counts_trim <- data.frame(t(data.frame(norm_counts)))                                                  
  norm_counts_trim[] <- lapply(norm_counts_trim, function(x) DescTools::Winsorize(x, probs = c(0, 0.97), type = 1))
  norm_counts_trim <- data.frame(t(norm_counts_trim))
  norm_counts_trim[norm_counts_trim == 0] <- 1
  norm_counts_trim <- log2(norm_counts_trim)
  colnames(norm_counts_trim) <- gsub("\\.", "-", colnames(norm_counts_trim))
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(norm_counts_trim, taxa_are_rows =  TRUE)
  return(ps_norm)
}

ps.16s_n <- wrench.norm.ps(ps.16s)
ps.its_n <- wrench.norm.ps(ps.its)

ps.16s.ctrl_n <- wrench.norm.ps(ps.16s.ctrl)
ps.its.ctrl_n <- wrench.norm.ps(ps.its.ctrl)
```

```R

```

```R

```

```R

```

```R

```

```R

```

```R

```

```R

```

```R

```

```R

```

```R

```
