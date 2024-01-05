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

