# Lettuce seedlings rapidly assemble their microbiome from the environment through deterministic processes

**Nesma Zakaria Mohamed, Leonardo Schena, Antonino Malacrinò**

## Abstract
Plant-associated microorganisms have a significant impact on plant biology, ecology, and evolution. Despite several studies that have examined the factors driving variation in plant microbiomes, the mechanisms behind the assembly of the plant microbiome are still poorly understood. In this study we use gnotobiotic plants to test (i) whether seedlings create a selective environment and drive the assembly of root and leaf microbiomes through deterministic or stochastic processes, and (ii) if seedlings structure the microbiome that is transferred through seeds using deterministic processes, and whether this pattern changes when seedlings are exposed to the environmental microbiome. Our results show that the microbiome of gnotobiotic plants (i.e., inherited through seeds) is not under the selective influence of the host plant, but this quickly changes when plants are exposed to soil microbiomes. Indeed, within one week, plants were able to select microorganisms from the inocula, assemble the root microbiome and from this, assemble the shoot microbiome. This study supports the hypothesis that plants at early developmental stages might exert strong selective activity on their microbiome, and contributes in clarifying the mechanisms of plant microbiome assembly.

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject numbers `PRJNA1029922` (16S amplicon metagenomics) and `PRJNA1029924` (ITS amplicon metagenomics).

Our pipeline included:
* DADA2 (Callahan et al. [2016](https://www.nature.com/articles/nmeth.3869))
* R  (R Core Team [2022](https://www.R-project.org/))

# Code

### **1.** [DADA2](/01_dada2)
Code for processing raw reads.

### **2.** [Data analysis](/02_data_analysis)
Data and code for reproducing all analyses.
