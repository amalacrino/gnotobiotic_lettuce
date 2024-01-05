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
library("ggridges")
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

## PERMANOVA - 16S

```R
sampledf <- data.frame(sample_data(ps.16s_n))
dist.mat <- phyloseq::distance(ps.16s_n, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ compartment, data = sampledf, permutations = perm)
pmv

RVAideMemoire::pairwise.perm.manova(dist.mat, sampledf$compartment, nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)
```

## NMDS - 16S

Figure 1, panels A and B

```R
plot.nmds <- function(ps, dist){
  dist.mat <- phyloseq::distance(ps, method = dist)
  cap_ord <- ordinate(physeq = ps, method = "NMDS", distance = dist.mat, formula = ~ 1)
  cap_plot <- plot_ordination(physeq = ps, ordination = cap_ord, axes = c(1,2)) +
    theme_bw(base_size = 20) +
    stat_ellipse(mapping = aes(fill = compartment),
                   alpha = 0.4,
                   geom = "polygon",
                   show.legend=T) +
    geom_point(mapping = aes(color = compartment), size = 5) +
    theme(legend.title= element_blank(), 
          legend.background = element_rect(color = NA),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          panel.grid = element_blank()) +
    scale_color_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot")) +
    scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))
}

plot1 <- plot.nmds(ps.16s_n, "unifrac")
plot2 <- plot.nmds(ps.16s_n, "wunifrac")

px <- ggarrange(plot1, plot2, ncol = 2,  align = "hv", common.legend = T)
px
```

## PERMANOVA - ITS

```R
sampledf <- data.frame(sample_data(ps.its_n))
dist.mat <- phyloseq::distance(ps.its_n, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ compartment, data = sampledf, permutations = perm)
pmv

RVAideMemoire::pairwise.perm.manova(dist.mat, sampledf$compartment, nperm = 999, progress = TRUE, p.method = "fdr", F = T, R2 = T)
```

## NMDS - ITS

Figure 1, panels C and D

```R
plot1 <- plot.nmds(ps.its_n, "unifrac")
plot2 <- plot.nmds(ps.its_n, "wunifrac")

px <- ggarrange(plot1, plot2, ncol = 2,  align = "hv", common.legend = T)
px
```

## Diversity - 16S

### Plot

Figure 2, panels A-D.

```R
div <- microbiome::alpha(ps.16s, index = c("observed", "diversity_shannon", "dominance_simpson"))
otus <- as.data.frame((otu_table(ps.16s)))
tree <- phy_tree(ps.16s)
div.pd <- pd(otus, tree, include.root = FALSE)
div.2 <- cbind(sample_data(ps.16s), div)
div.2 <- cbind(div.2, div.pd)

div_plot1 <- ggplot(div.2, aes(x = compartment, y = PD, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Phylogenetic diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 15) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot2 <- ggplot(div.2, aes(x = compartment, y = dominance_simpson, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Simpson's dominance") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 0.05) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot3 <- ggplot(div.2, aes(x = compartment, y = diversity_shannon, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Shannon's diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 7) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot4 <- ggplot(div.2, aes(x = compartment, y = observed, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Observed richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 400) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

px <- ggarrange(div_plot1, div_plot3, div_plot2, div_plot4,
                ncol = 4,  align = "hv")
px
```

### Test

Figure 2 and Table S1.

#### Model PD
```R
model <- lm(PD ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Shannon
```R
model <- lm(diversity_shannon ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Simpson
```R
model <- lm(dominance_simpson ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Observed
```R
model <- lm(observed ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

## Diversity - ITS

### Plot

Figure 2, panels E-H.

```R
div <- microbiome::alpha(ps.its, index = c("observed", "diversity_shannon", "dominance_simpson"))
otus <- as.data.frame((otu_table(ps.its)))
tree <- phy_tree(ps.its)
div.pd <- pd(otus, tree, include.root = FALSE)
div.2 <- cbind(sample_data(ps.its), div)
div.2 <- cbind(div.2, div.pd)

div_plot1 <- ggplot(div.2, aes(x = compartment, y = PD, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Phylogenetic diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 40) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot2 <- ggplot(div.2, aes(x = compartment, y = dominance_simpson, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Simpson's dominance") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 1) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot3 <- ggplot(div.2, aes(x = compartment, y = diversity_shannon, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Shannon's diversity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 6) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))

div_plot4 <- ggplot(div.2, aes(x = compartment, y = observed, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  labs(y = "Observed richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0, 300) +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot"))


px <- ggarrange(div_plot1, div_plot3, div_plot2, div_plot4,
                ncol = 4,  align = "hv")
px
```

### Test

Figure 2 and Table S1.

#### Model PD
```R
model <- lm(PD ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Shannon
```R
model <- lm(diversity_shannon ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Simpson
```R
model <- lm(dominance_simpson ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

#### Model Observed
```R
model <- lm(observed ~ compartment, data = div.2)
Anova(model)
m1 <- emmeans(model, "compartment")
pairs(m1)
multcomp::cld(object = m1, Letters = letters)
```

## Differential abundance of microbial taxa - 16S

Figure 3, panels A-C.

```R
diff.taxa.fun <- function(ps, g1, g2){
  ks <- sample_data(ps)[["compartment"]] %in% c(g1, g2)
  ps_da <- prune_samples(samples = ks, ps)
  mas_1 <- Maaslin2(
    input_data = data.frame(t(otu_table(ps_da))),
    input_metadata = data.frame(sample_data(ps_da)),
    output = "demo_output6", 
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.05,
    fixed_effects = "compartment",
    random_effects = c("group"),
    correction = "BH",
    standardize = FALSE,
    cores = 4)
  df.diff <- mas_1$results
  tax.table <- as.data.frame(tax_table(ps))
  tax.table <- setDT(tax.table, keep.rownames = TRUE)[]
  tx <- merge(df.diff, tax.table, by.x = "feature", by.y = "rn")
  return(tx)
}

plot.diff.taxa <- function(df, g1, g2, title){
  df.diff <- df
  df.diff$diffexpressed <- "no changes"
  df.diff$diffexpressed[df.diff$coef > 0 & df.diff$qval < 0.05] <- paste0(g2)
  df.diff$diffexpressed[df.diff$coef < 0 & df.diff$qval < 0.05] <- paste0(g1)
  plot <- ggplot(data=df.diff) +
        theme_bw(base_size = 12) +
        geom_point(aes(x = coef, y = -log10(qval), colour = diffexpressed)) +
        scale_color_manual(values=c("blue", "red", "black")) +
        geom_vline(xintercept=0, col="black", linetype = "longdash") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "longdash") +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", 
                                              colour = "white", 
                                              size = 0.5, linetype = "solid"),
              legend.position="none",
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        ggtitle(paste0(title)) +
        xlab(expression(paste(Log[2], " Fold Changes"))) +
        ylab(expression(paste(-Log[10], " P"))) +
        scale_color_manual(name = "Legend", values=c("#4daf4a", "#e41a1c", "#000000"), labels = c(g2, g1), breaks =c(g2, g1))
  return(plot)
}

df.da.IR.16s <- diff.taxa.fun(ps.16s_n, "inocula", "root")
df.da.IS.16s <- diff.taxa.fun(ps.16s_n, "inocula", "shoot")
df.da.RS.16s <- diff.taxa.fun(ps.16s_n, "root", "shoot")

plot1 <- plot.diff.taxa(df.da.IR.16s, "inocula", "root", "root vs inocula")
plot2 <- plot.diff.taxa(df.da.IS.16s, "inocula", "shoot", "shoot vs inocula")
plot3 <- plot.diff.taxa(df.da.RS.16s, "root", "shoot", "shoot vs root")

px <- ggarrange(plot1, plot2, plot3, ncol = 3,  align = "hv")
px
```

## Differential abundance of microbial taxa - ITS

Figure 3, panels D-F.

```R
df.da.IR.its <- diff.taxa.fun(ps.its_n, "inocula", "root")
df.da.IS.its <- diff.taxa.fun(ps.its_n, "inocula", "shoot")
df.da.RS.its <- diff.taxa.fun(ps.its_n, "root", "shoot")

plot1 <- plot.diff.taxa(df.da.IR.its, "inocula", "root", "root vs inocula")
plot2 <- plot.diff.taxa(df.da.IS.its, "inocula", "shoot", "shoot vs inocula")
plot3 <- plot.diff.taxa(df.da.RS.its, "root", "shoot", "shoot vs root")

px <- ggarrange(plot1, plot2, plot3, ncol = 3,  align = "hv")
px
```

## ASVs shared between compartments - 16S

Figure 4, panel A.

### Plot

```R
ps_test <- ps.16s

get.numbers <- function(ps, group){
  temp <- sample_data(ps)[["group"]] %in% group
  ps.filt <- prune_samples(samples = temp, ps)
  '%ni%' <- Negate('%in%')
  
  temp <- sample_data(ps.filt)[["compartment"]] %in% "root"
  ps.filt.r <- prune_samples(samples = temp, ps.filt)
  temp <- sample_data(ps.filt)[["compartment"]] %in% "shoot"
  ps.filt.s <- prune_samples(samples = temp, ps.filt)
  temp <- sample_data(ps.filt)[["compartment"]] %in% "inocula"
  ps.filt.i <- prune_samples(samples = temp, ps.filt)
  
  ps.filt.r <- filter_taxa(ps.filt.r, function (x) {sum(x > 0) > 0}, prune=TRUE)
  ps.filt.s <- filter_taxa(ps.filt.s, function (x) {sum(x > 0) > 0}, prune=TRUE)
  ps.filt.i <- filter_taxa(ps.filt.i, function (x) {sum(x > 0) > 0}, prune=TRUE)
  
  ps.filt.r <- t(as.data.frame(otu_table(ps.filt.r)))
  ps.filt.s <- t(as.data.frame(otu_table(ps.filt.s)))
  ps.filt.i <- t(as.data.frame(otu_table(ps.filt.i)))
  
  ps.filt.r <- row.names(ps.filt.r)
  ps.filt.s <- row.names(ps.filt.s)
  ps.filt.i <- row.names(ps.filt.i)
  

  g <- Reduce(intersect, list(ps.filt.i, ps.filt.r, ps.filt.s))
  a <- ps.filt.i %ni% unique(c(ps.filt.r, ps.filt.s))
  b <- ps.filt.r %ni% unique(c(ps.filt.i, ps.filt.s))
  c <- ps.filt.s %ni% unique(c(ps.filt.i, ps.filt.r))
  
  d <- Reduce(intersect, list(ps.filt.r, ps.filt.s)) %ni% g
  e <- Reduce(intersect, list(ps.filt.s, ps.filt.i)) %ni% g
  f <- Reduce(intersect, list(ps.filt.r, ps.filt.i)) %ni% g
  

  res <- data.frame("inoculum" = length(a[a == TRUE]),
                    "roots" = length(b[b == TRUE]),
                    "shoot"= length(c[c == TRUE]),
                    "roots_shoot" = length(d[d == TRUE]),
                    "shoot_inoculum"= length(e[e == TRUE]),
                    "roots_inoculum"  = length(f[f == TRUE]),
                    "all_compartments" = length(g))
  return(res)
}

sampledf <- data.frame(sample_data(ps_test))
list.inocula <- unique(sampledf$group)
model_calculator <- sapply(list.inocula, get.numbers, ps = ps_test, simplify = FALSE, USE.NAMES = TRUE)
res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
colnames(res)[1] <- "sample_id"
library(tidyverse)
new_data <- melt(res, id='sample_id')
new_dat1 <- new_data %>% 
  arrange(sample_id)
colnames(new_dat1)[3] <- "observed"

generate.rand.ps <- function(psobject){
  psx <- psobject
  ps <- as.data.frame(otu_table(psx))
  set.seed(100)
  otu <- randomizeMatrix(ps,null.model = "frequency",iterations = 1000)
  otu_table(psx) <- otu_table(otu, taxa_are_rows = F)
  return(psx)
}

ps.random <- generate.rand.ps(ps_test)
sampledf <- data.frame(sample_data(ps.random))
list.inocula <- unique(sampledf$group)
model_calculator <- sapply(list.inocula, get.numbers, ps = ps.random, simplify = FALSE, USE.NAMES = TRUE)
res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
colnames(res)[1] <- "sample_id"

new_data <- melt(res, id='sample_id')
new_dat2 <- new_data %>% 
  arrange(sample_id)
colnames(new_dat2)[3] <- "random"

nd_df <- new_dat1
nd_df$random <- new_dat2$random

density_plot1 <- ggplot(nd_df, aes(x = observed, y = variable)) +
  theme_bw(base_size = 15) +
  geom_density_ridges(fill = "#a1d76a", rel_min_height = 0.01, alpha = 0.7, scale = 0.8,
                      jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  geom_density_ridges(aes(x = random, y = variable), fill = "#e9a3c9", rel_min_height = 0.01, alpha = 0.7, scale = 0.8,
                      jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '*', point_size = 3, point_alpha = 1, alpha = 0.7) +
  labs(x = "# of ASVs",
       y = "intersections") +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  xlim(0, 300)
density_plot1
```

### Test

```R
list.mol <- unique(nd_df$variable)
model_calculator.1 <- sapply(list.mol,  
                             function(x){
                               df <- nd_df %>% filter(variable == x) %>% 
                                 melt(id.vars = c("sample_id", "variable")) %>% 
                                 dplyr::rename("compartment" = "variable", "group" = "variable.1")
                               model <- glmer(value ~ group * (1|sample_id), data = df)
                               aaa <-  Anova(model)
                               aaa$sig = c(rep('',length(aaa$`Pr(>Chisq)`)))
                               makeStars <- function(x){
                                 stars <- c("****", "***", "**", "*", "ns")
                                 vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
                                 i <- findInterval(x, vec)
                                 stars[i] }
                               aaa$sig <- makeStars(aaa$`Pr(>Chisq)`)
                               aaa <- aaa[1,]
                               row.names(aaa) <- x
                               return(aaa)},
                             simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator.1)
res <- setDT(res, keep.rownames = TRUE)[]
res
```

## ASVs shared between compartments - ITS

Figure 4, panel B.

### Plot

```R
ps_test <- ps.its
sampledf <- data.frame(sample_data(ps_test))
sel <- sampledf %>% group_by(group) %>% summarise(n = n()) %>% filter(n > 2)
temp <- sample_data(ps_test)[["group"]] %in% sel$group
ps_test <- prune_samples(samples = temp, ps_test)
sampledf <- data.frame(sample_data(ps_test))

list.inocula <- unique(sampledf$group)
model_calculator <- sapply(list.inocula, get.numbers, ps = ps_test, simplify = FALSE, USE.NAMES = TRUE)
res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
colnames(res)[1] <- "sample_id"

new_data <- melt(res, id='sample_id')
new_dat1 <- new_data %>% 
  arrange(sample_id)
colnames(new_dat1)[3] <- "observed"


ps.random <- generate.rand.ps(ps_test)
sampledf <- data.frame(sample_data(ps.random))
list.inocula <- unique(sampledf$group)
model_calculator <- sapply(list.inocula, get.numbers, ps = ps.random, simplify = FALSE, USE.NAMES = TRUE)
res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
colnames(res)[1] <- "sample_id"

new_data <- melt(res, id='sample_id')
new_dat2 <- new_data %>% 
  arrange(sample_id)
colnames(new_dat2)[3] <- "random"

nd_df <- new_dat1
nd_df$random <- new_dat2$random

density_plot2 <- ggplot(nd_df, aes(x = observed, y = variable)) +
  theme_bw(base_size = 15) +
  geom_density_ridges(fill = "#a1d76a", rel_min_height = 0.01, alpha = 0.7, scale = 0.8,
                      jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  geom_density_ridges(aes(x = random, y = variable), fill = "#e9a3c9", rel_min_height = 0.01, alpha = 0.7, scale = 0.8,
                      jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0), point_shape = '*', point_size = 3, point_alpha = 1, alpha = 0.7) +
  labs(x = "# of ASVs",
       y = "intersections") +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  xlim(0, 300)
density_plot2
```

### Test

```R
model_calculator.1 <- sapply(list.mol,  
                             function(x){
                               df <- nd_df %>% filter(variable == x) %>% 
                                 melt(id.vars = c("sample_id", "variable")) %>% 
                                 rename("compartment" = "variable", "group" = "variable.1")
                               model <- glmer(value ~ group * (1|sample_id), data = df)
                               aaa <-  Anova(model)
                               aaa$sig = c(rep('',length(aaa$`Pr(>Chisq)`)))
                               makeStars <- function(x){
                                 stars <- c("****", "***", "**", "*", "ns")
                                 vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
                                 i <- findInterval(x, vec)
                                 stars[i] }
                               aaa$sig <- makeStars(aaa$`Pr(>Chisq)`)
                               aaa <- aaa[1,]
                               row.names(aaa) <- x
                               return(aaa)},
                             simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator.1)
res <- setDT(res, keep.rownames = TRUE)[]
res
```

## beta-NTI - 16S

### Plot and test

Figure 5, panel A.

```R
comm <- as.data.frame(t(otu_table(ps.16s.ctrl_n)))
phy <- phy_tree(ps.16s.ctrl_n)
phy.dist <- cophenetic(phy)
comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
ks <- sample_data(ps.16s.ctrl_n)
bnti <- cbind(ks, comm.sesmpd)

betapart_plot1 <- ggplot(bnti, aes(x = compartment, y = -mntd.obs.z, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red") +
  labs(y = "beta-NTI") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root gnotobiotic", "shoot gnotobiotic"), breaks = c("inocula", "root gnotobiotic", "shoot gnotobiotic")) +
  ylim(-2.5,9)
betapart_plot1
```

Figure 5, panel C.

```R
comm <- as.data.frame(t(otu_table(ps.16s_n)))
phy <- phy_tree(ps.16s_n)
phy.dist <- cophenetic(phy)
comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
ks <- sample_data(ps.16s_n)
bnti <- cbind(ks, comm.sesmpd)
bnti <- bnti[which(bnti$compartment != "inocula"),]

betapart_plot1 <- ggplot(bnti, aes(x = compartment, y = -mntd.obs.z, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red") +
  labs(y = "beta-NTI") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot")) +
  ylim(-2.5, 9)
betapart_plot1
```

## beta-NTI - ITS

### Plot

Figure 5, panel B.

```R
comm <- as.data.frame(t(otu_table(ps.its.ctrl_n)))
phy <- phy_tree(ps.its.ctrl_n)
phy.dist <- cophenetic(phy)
comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
ks <- sample_data(ps.its.ctrl_n)
bnti <- cbind(ks, comm.sesmpd)

betapart_plot1 <- ggplot(bnti, aes(x = compartment, y = -mntd.obs.z, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red") +
  labs(y = "beta-NTI") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root gnotobiotic", "shoot gnotobiotic"), breaks = c("inocula", "root gnotobiotic", "shoot gnotobiotic")) +
  ylim(-2.5,9)
betapart_plot1
```

Figure 5, panel D.

```R
comm <- as.data.frame(t(otu_table(ps.its_n)))
phy <- phy_tree(ps.its_n)
phy.dist <- cophenetic(phy)
comm.sesmpd <- ses.mntd(comm, phy.dist, null.model = "richness", abundance.weighted = T, runs = 999)
ks <- sample_data(ps.its_n)
bnti <- cbind(ks, comm.sesmpd)
bnti <- bnti[which(bnti$compartment != "inocula"),]

betapart_plot1 <- ggplot(bnti, aes(x = compartment, y = -mntd.obs.z, fill = compartment)) +
  theme_bw(base_size = 15) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black") +
  geom_hline(yintercept=2, linetype="dashed", color = "red") +
  geom_hline(yintercept=-2, linetype="dashed", color = "red") +
  labs(y = "beta-NTI") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", face = "italic"),
        legend.position = "none") +
  scale_fill_manual(name = "Legend", values=c("#7570b3", "#d95f02", "#1b9e77"), labels = c("inocula", "root", "shoot"), breaks = c("inocula", "root", "shoot")) +
  ylim(-2.5, 9)
betapart_plot1
```

## Taxa plots

### 16S

Figure S2.

```R
glom <- microbiome::aggregate_taxa(ps.16s_n, "Genus")
glom <- microbiome::transform(glom, "compositional")
dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarise(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(OTU %in% filt.gen$Genus))
dat <- dat %>% group_by(Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 


dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarise(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(OTU %in% filt.gen$Genus))
dat <- dat %>% group_by(compartment, Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 

nb.cols <- length(unique(dat$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

taxa_plot <- ggplot(dat, aes(x = as.factor(compartment), y = cs, fill = Genus)) +
                        theme_bw(base_size = 14) +
                        geom_bar(stat="identity") +
                        labs(y = "Relative proportion") +
                        theme(legend.background = element_rect(fill="white"),
                              legend.key = element_rect(fill="transparent"),
                              legend.text = element_text(size = 12),
                              axis.text.x = element_text(color="black"),
                              axis.text.y = element_text(color="black"),
                              panel.grid = element_blank()) +
                        scale_y_continuous(labels = scales::percent) +
                        scale_fill_manual(values = mycolors) +
                        labs(y = "Relative abundance", x="")
taxa_plot
```

### ITS

Figure S3.

```R
glom <- microbiome::aggregate_taxa(ps.its_n, "Genus")
glom <- microbiome::transform(glom, "compositional")
dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarise(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(OTU %in% filt.gen$Genus))
dat <- dat %>% group_by(Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 


dat <- psmelt(glom)
filt.gen <- dat %>% group_by(Genus) %>% summarise(mean = mean(Abundance)) %>% filter(mean <= 0.01)
dat <- subset(dat, !(OTU %in% filt.gen$Genus))
dat <- dat %>% group_by(compartment, Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 

nb.cols <- length(unique(dat$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

taxa_plot <- ggplot(dat, aes(x = as.factor(compartment), y = cs, fill = Genus)) +
                        theme_bw(base_size = 14) +
                        geom_bar(stat="identity") +
                        labs(y = "Relative proportion") +
                        theme(legend.background = element_rect(fill="white"),
                              legend.key = element_rect(fill="transparent"),
                              legend.text = element_text(size = 12),
                              axis.text.x = element_text(color="black"),
                              axis.text.y = element_text(color="black"),
                              panel.grid = element_blank()) +
                        scale_y_continuous(labels = scales::percent) +
                        scale_fill_manual(values = mycolors) +
                        labs(y = "Relative abundance", x="") 
ggsave(taxa_plot, filename = "figures/taxa_its.pdf", dpi = 600,  width = 5, height = 6, units = "in")
taxa_plot
```

### Gnotobiotic plants

Figure S4, panel A.

```R
glom <- microbiome::transform(ps.16s.ctrl_n, "compositional")
dat.ctrl <- psmelt(glom)
dat.ctrl1 <- dat.ctrl %>% group_by(Sample, Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 

nb.cols <- length(unique(dat.ctrl1$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
taxa_plot <- ggplot(dat.ctrl1, aes(x = as.factor(Sample), y = cs, fill = Genus)) +
                        theme_bw(base_size = 14) +
                        geom_bar(stat="identity") +
                        labs(y = "Relative proportion") +
                        theme(legend.background = element_rect(fill="white"),
                              legend.key = element_rect(fill="transparent"),
                              legend.text = element_text(size = 12),
                              axis.text.x = element_text(color="black"),
                              axis.text.y = element_text(color="black"),
                              panel.grid = element_blank()) +
                        scale_y_continuous(labels = scales::percent) +
                        scale_fill_manual(values = mycolors) +
                        labs(y = "Relative abundance", x="") 
taxa_plot
```

Figure S4, panel B.

```R
glom <- microbiome::transform(ps.its.ctrl_n, "compositional")
dat.ctrl <- psmelt(glom)
dat.ctrl1 <- dat.ctrl %>% group_by(Sample, Genus) %>% summarise(cs = mean(Abundance)) %>% mutate(cs = cs/sum(cs)) 

nb.cols <- length(unique(dat.ctrl1$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
taxa_plot <- ggplot(dat.ctrl1, aes(x = as.factor(Sample), y = cs, fill = Genus)) +
                        theme_bw(base_size = 14) +
                        geom_bar(stat="identity") +
                        labs(y = "Relative proportion") +
                        theme(legend.background = element_rect(fill="white"),
                              legend.key = element_rect(fill="transparent"),
                              legend.text = element_text(size = 12),
                              axis.text.x = element_text(color="black"),
                              axis.text.y = element_text(color="black"),
                              panel.grid = element_blank()) +
                        scale_y_continuous(labels = scales::percent) +
                        scale_fill_manual(values = mycolors) +
                        labs(y = "Relative abundance", x="")
taxa_plot
```
