setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
setwd("C:/Users/jas/Google Drive/Metagenomics/oab/")
setwd("/Users/administrator/Documents/jaspreet/oab/oab")

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

# Make phyloseq object ----------------------------------------------------
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

#phyloseq object

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

d_r = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                    Phylum == "p:Actinobacteria")

###############Soil------------------------------------

###make pc soil phyloseq

decon = subset_samples(d, Source == "soil")
decon
decon = subset_samples(decon, Species == "P. cooperi")
decon
decon = subset_samples(decon, Pop_size != "Con")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency
source("scripts/decontaminate_phyloseq.R")

d.pc = subset_samples(decon.d, Sample_or_Control == "Sample")
d.pc
d.pc = prune_taxa(taxa_sums(d.pc) >= 1, d.pc)
d.pc
comb = tax_glom(d.pc, "Kingdom")
comb
taxa_sums(comb)

###make pp soil phyloseq

d.pp = subset_samples(d, Source == "soil")
d.pp
d.pp = subset_samples(d.pp, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp

####Combine phyloseq objects

d.fin = merge_phyloseq(d_r, d.pc, d.pp)
d.fin
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source,
                               ".", sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

a = merge_samples(d.fin, "int")
a

b = subset_taxa(a, Family == "f:Burkholderiaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Pseudomonadaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Bacillales_incertae_sedis")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Halomonadaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Xanthomonadaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Micrococcaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

