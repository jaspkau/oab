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

###Filter d.r.net for most abundant bacterial phyla observed in roots

d.fin = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                      Phylum == "p:Actinobacteria")

d.fin = subset_taxa(d.fin, Family == "f:Pseudomonadaceae")

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$replicate)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

# NETWORK ANALYSIS -------------------------------------------------------
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)

##make root otu file
##"Remove taxa not seen more than 3 times in at least 5% of the samples. 
#This protects against an OTU with small mean & trivially large C.V.
#d.r.net = filter_taxa(d.r.net, function(x) sum(x > 1) > (0.05*length(x)), TRUE)

d.r.net = subset_samples(d.fin, Species == "P. cooperi")
d.r.net
d.r.net = merge_samples(d.r.net, "int")
otu.r.net = data.frame(otu_table(d.r.net)) ##it should be non-normalized
row.names(otu.r.net) = seq(1:6)
colnames(otu.r.net) = paste(gsub("otu", "pc", colnames(otu.r.net)))
taxa_names(d.r.net) = paste(gsub("otu", "pc", taxa_names(d.r.net)))

##make soil otu file
d.s.net = subset_samples(d.fin, Species == "P. praeclara")
d.s.net = merge_samples(d.s.net, "int")
d.s.net
otu.s.net = data.frame(otu_table(d.s.net)) ##it should no non-normalized
row.names(otu.s.net) = seq(1:8)
colnames(otu.s.net) = paste(gsub("otu", "pp", colnames(otu.s.net)))
taxa_names(d.s.net) = paste(gsub("otu", "pp", taxa_names(d.s.net)))

d.net = merge_phyloseq(d.r.net, d.s.net)

###Sparcc
##merge root and soil otu files together

net = merge(otu.s.net, otu.r.net, by = "row.names", all = TRUE)
net = rbind(otu.r.net, otu.s.net)
net = t(net)
write.csv(net, file = "results/net.csv")
#convert csv to .txt file
#and then do analyses with spaccWrapper.sh script and then import the files with following codes for editing to use in Cytoscape

library(reshape2)
spec.cor = read.delim("results/net_pseudo_species/sim_cor.txt", sep = "\t")
row.names(spec.cor) = spec.cor[,1]
spec.cor = melt(spec.cor)

pval = read.delim("results/net_pseudo_species/pvals_two_sided.txt", sep = "\t")
row.names(pval) = pval[,1]
pval = melt(pval)
library(splitstackshape)
library(plyr)
pval = plyr::rename(pval, c("value" = "p"))
spec.cor$p = pval$p

spec.cor.t.sel = spec.cor[grepl("pp", spec.cor$X), ] #target selection
spec.cor.fin.sel = spec.cor.t.sel[grepl("pc", spec.cor.t.sel$variable), ]
write.csv(spec.cor.fin.sel, file = "results/net_pseudo_species/spec.cor.csv")

###select for co-abundance network
spec.cor.fin.sel1 = spec.cor.fin.sel[spec.cor.fin.sel$value > 0.6,]
spec.cor.fin.sel1 = spec.cor.fin.sel1[spec.cor.fin.sel1$p == 0,]
write.csv(spec.cor.fin.sel1, file = "results/net_pseudo/spec.cor.fin.sel.csv")

###select for co-exclusion network
spec.cor.fin.sel2 = spec.cor.fin.sel[spec.cor.fin.sel$value < -0.6,]
spec.cor.fin.sel2 = spec.cor.fin.sel2[spec.cor.fin.sel2$p == 0,]
write.csv(spec.cor.fin.sel2, file = "results/net_pseudo/spec.cor.fin.sel_co_exclu.csv")

###select pair of root and soil OTUs
net = read.delim("results/net_pseudo_species/spec.cor.csv", sep = ",")
net$X = gsub("pp", "", net$X)
net$variable = gsub("pc", "", net$variable)
net$sel = ifelse((net$X == net$variable)==TRUE, paste("sel"), paste(""))
#or
net$sel = ifelse(net$X == net$variable, paste("sel"), paste(""))
net = subset(net, sel == "sel")
write.csv(net, file = "results/net_pseudo_species/net_sel.csv")
