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

comb = tax_glom(d_r, "Kingdom")
comb
taxa_sums(comb)

####MOST ABUDANT PHYLA IN ROOTS

a = subset_taxa(d_r, Phylum == "p:Proteobacteria")
a
b = tax_glom(a, "Phylum")
b
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d_r, Phylum == "p:Firmicutes")
a
b = tax_glom(a, "Phylum")
b
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d_r, Phylum == "p:Actinobacteria")
a
b = tax_glom(a, "Phylum")
b
x = taxa_sums(b)/taxa_sums(comb)
x

d_r = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                      Phylum == "p:Actinobacteria")

####OTU shared between species and unique ones

d.pc = subset_samples(d_r, Species == "P. cooperi")
d.pc
d.pc = prune_taxa(taxa_sums(d.pc) >= 1, d.pc)
d.pc
d.pp = subset_samples(d_r, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp
d.shared = prune_taxa(taxa_names(d.pp) %in% taxa_names(d.pc), d.pp)
d.shared

##############phenological stages
sample_data(d_r)$int = paste(sample_data(d_r)$Species,".",sample_data(d_r)$Stage)
sample_data(d_r)$int = gsub(" ", "", sample_data(d_r)$int)

a = merge_samples(d_r, "int")
a

comb = tax_glom(a, "Kingdom")
comb
sample_sums(comb)

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

b = subset_taxa(a, Family == "f:Enterobacteriaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Xanthomonadaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Bacillales_incertae_sedis")
b
b = tax_glom(b, "Family")
b
sample_sums(b)


###############Soil--------------------------------------------------------

###make pc soil phyloseq

decon = subset_samples(d, Source == "soil")
decon
decon = subset_samples(decon, Species == "P. cooperi")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency
source("scripts/decontaminate_phyloseq.R")

d.pc = subset_samples(decon.d, Sample_or_Control == "Sample")
d.pc
d.pc = prune_taxa(taxa_sums(d.pc) >= 1, d.pc)
d.pc
comb.pc = tax_glom(d.pc, "Kingdom")
comb.pc
taxa_sums(comb.pc)

###make pp soil phyloseq

d.pp = subset_samples(d, Source == "soil")
d.pp
d.pp = subset_samples(d.pp, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp

comb.pp = tax_glom(d.pp, "Kingdom")
comb.pp
taxa_sums(comb.pp)

####Combine phyloseq objects

d.s = merge_phyloseq(d.pc, d.pp)
d.s

comb = tax_glom(d.s, "Kingdom")
comb
taxa_sums(comb)

a = subset_taxa(d.s, Phylum == "p:Proteobacteria")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d.s, Phylum == "p:Acidobacteria")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d.s, Phylum == "p:Bacteroidetes")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d.s, Phylum == "p:Verrucomicrobia")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d.s, Phylum == "p:Planctomycetes")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

a = subset_taxa(d.s, Phylum == "p:Actinobacteria")
a
b = tax_glom(a, "Phylum")
b
taxa_sums(b)
x = taxa_sums(b)/taxa_sums(comb)
x

#####after selecting for bacterial phyla -------------------------

d.fin = subset_taxa(d.s, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")

d.pc3 = subset_samples(d.fin, Species == "P. cooperi")
d.pc3
d.pc3 = prune_taxa(taxa_sums(d.pc3) >= 1, d.pc3)
d.pc3

d.pp3 = subset_samples(d.fin, Species == "P. praeclara")
d.pp3
d.pp3 = prune_taxa(taxa_sums(d.pp3) >= 1, d.pp3)
d.pp3

d.pc.s_pp.s = prune_taxa(taxa_names(d.pp3) %in% taxa_names(d.pc3), d.pp3)
d.pc.s_pp.s

##########root OTUs in soil#######################
##########################################

d.r_d.s = prune_taxa(taxa_names(d.fin) %in% taxa_names(d_r), d.fin)
d.r_d.s
d.r_pcs = prune_taxa(taxa_names(d.pc3) %in% taxa_names(d_r), d.pc3)
d.r_pcs
comb = tax_glom(d.r_pcs, "Kingdom")
comb
taxa_sums(comb)/taxa_sums(comb.pc)

d.r_pps = prune_taxa(taxa_names(d.pp3) %in% taxa_names(d_r), d.pp3)
d.r_pps
comb = tax_glom(d.r_pps, "Kingdom")
comb
taxa_sums(comb)
taxa_sums(comb)/taxa_sums(comb.pp)

####################################Pop_size counts
####Combine phyloseq objects

d.fin = merge_phyloseq(d_r, d.pc, d.pp)
d.fin
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source,
                               ".", sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source)

a = merge_samples(d.fin, "int")
a

comb = tax_glom(a, "Kingdom")
comb
sample_sums(comb)

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

b = subset_taxa(a, Family == "f:Shewanellaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

b = subset_taxa(a, Family == "f:Rhizobiaceae")
b
b = tax_glom(b, "Family")
b
sample_sums(b)

