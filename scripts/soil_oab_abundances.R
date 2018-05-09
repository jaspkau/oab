setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
setwd("C:/Users/jaspr/Google Drive/Metagenomics/oab/")
setwd("/Users/administrator/Documents/jaspreet/oab/oab/")

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

###ROOT OMF ANALYSIS......................................
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "root")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence
source("scripts/decontaminate_phyloseq.R")

d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

r.otus = taxa_names(d_r)

###make soil phyloseq

decon = subset_samples(d, Source == "soil")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence
source("scripts/decontaminate_phyloseq.R")

d_s = subset_samples(decon.d, Sample_or_Control == "Sample")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

otu_s = data.frame(otu_table(d_s))
otu_s2 = otu_s[r.otus,]
otu_s2 = otu_s2[,colSums(otu_s2) > 0]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_oab = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_oab
d_oab = prune_taxa(taxa_sums(d_oab) >= 1, d_s)
d_oab
