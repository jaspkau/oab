setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")
setwd("/Users/administrator/Desktop/pico/")
setwd("/Users/administrator/Documents/jaspreet/oab/oab/")

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

library("Biostrings")
library("ggplot2")
#devtools::install_github("GuangchuangYu/treeio")
#devtools::install_github("GuangchuangYu/ggtree")
library(treeio)
library(ape)
library(ggtree)

# Make phyloseq object ----------------------------------------------------
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "root")
decon
#decon = subset_samples(decon, Population2 != "X")
#decon

####decontaminate phyloseq object based on frequency
source("scripts/decontaminate_phyloseq.R")

d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

d.fin = subset_taxa(d_r, Family == "f:Pseudomonadaceae")
taxa_names(d.fin) = gsub("otu", "denovo", taxa_names(d.fin))

#####
#Bring in the tree
#####
source("http://www.phytools.org/read.newick/v0.5/read.newick.R")
tre = read.newick(file = "results/raxml_pseudo/Untitled")

d.test = merge_phyloseq(d.fin, tre)
d.test = prune_samples(sample_sums(d.test) > 0, d.test)
unif_w = UniFrac(d.test, weighted=T, normalized = T)
unif_w = as.matrix(unif_w)
wuf = as.dist(unif_w)
#wuf[is.na(wuf)] <- 1

a = adonis(wuf ~ sample_data(d.test)$Species + sample_data(d.test)$Stage)
