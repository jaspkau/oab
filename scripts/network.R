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

###Filter d.fin for most abundant bacterial phyla observed in roots

d.fin = subset_taxa(d.fin, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")

d.fin = subset_taxa(d.fin, Family == "f:Pseudomonadaceae")

# NETWORK ANALYSIS -------------------------------------------------------
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)

##make root otu file
##"Remove taxa not seen more than 3 times in at least 5% of the samples. 
#This protects against an OTU with small mean & trivially large C.V.
#d.r.net = filter_taxa(d.fin, function(x) sum(x > 1) > (0.05*length(x)), TRUE)

d.r.net = subset_samples(d.fin, Source == "root")
d.r.net
d.r.net = prune_taxa(taxa_sums(d.r.net) >= 1, d.r.net)
d.r.net
d.r.net = merge_samples(d.r.net, "Population")
otu.r.net = data.frame(otu_table(d.r.net)) ##it should be non-normalized
###slelect OTUs first from which you want to build network, rows shud be samples
#otu.r.net = otu.r.net[row.names(otu.r.net) %in% row.names(sim.kw.popsize),]
colnames(otu.r.net) = paste(gsub("otu", "r", colnames(otu.r.net)))
taxa_names(d.r.net) = paste(gsub("otu", "r", taxa_names(d.r.net)))

d.s = subset_samples(d.fin, Pop_size == "L"|Pop_size == "S")

##make soil otu file
d.s.net = subset_samples(d.s, Source == "soil")
d.s.net = prune_taxa(taxa_sums(d.s.net) >= 1, d.s.net)
d.s.net
#d.s.net = filter_taxa(d.s.net, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
d.s.net = merge_samples(d.s.net, "Population")
d.s.net
otu.s.net = data.frame(otu_table(d.s.net)) ##it should no non-normalized
colnames(otu.s.net) = paste(gsub("otu", "s", colnames(otu.s.net)))
taxa_names(d.s.net) = paste(gsub("otu", "s", taxa_names(d.s.net)))

d.net = merge_phyloseq(d.r.net, d.s.net)

###Sparcc
##merge root and soil otu files together

net = cbind(otu.s.net, otu.r.net, by = "row.names")
net = net[,-ncol(net)]

net = t(net)
write.csv(net, file = "results/net.csv")
#convert csv to .txt file
#and then do analyses with spaccWrapper.sh script and then import the files with following codes for editing to use in Cytoscape

library(reshape2)
spec.cor = read.delim("results/net_pseudo/sim_cor.txt", sep = "\t")
row.names(spec.cor) = spec.cor[,1]
spec.cor = melt(spec.cor)

pval = read.delim("results/net_pseudo/pvals_two_sided.txt", sep = "\t")
row.names(pval) = pval[,1]
pval = melt(pval)
library(splitstackshape)
library(plyr)
pval = plyr::rename(pval, c("value" = "p"))
spec.cor$p = pval$p

spec.cor.t.sel = spec.cor[grepl("s", spec.cor$X), ] #target selection
spec.cor.fin.sel = spec.cor.t.sel[grepl("r", spec.cor.t.sel$variable), ]
write.csv(spec.cor.fin.sel, file = "results/net_pseudo/spec.cor.csv")

###select for co-abundance network
spec.cor.fin.sel1 = spec.cor.fin.sel[spec.cor.fin.sel$value > 0.6,]
spec.cor.fin.sel1 = spec.cor.fin.sel1[spec.cor.fin.sel1$p == 0,]
write.csv(spec.cor.fin.sel1, file = "results/net_pseudo/spec.cor.fin.sel.csv")

###select for co-exclusion network
spec.cor.fin.sel2 = spec.cor.fin.sel[spec.cor.fin.sel$value < -0.6,]
spec.cor.fin.sel2 = spec.cor.fin.sel2[spec.cor.fin.sel2$p == 0,]
write.csv(spec.cor.fin.sel2, file = "results/net_pseudo/spec.cor.fin.sel_co_exclu.csv")

###select pair of root and soil OTUs
net = read.delim("results/net_pseudo/spec.cor.csv", sep = ",")
net$X = gsub("s", "", net$X)
net$variable = gsub("r", "", net$variable)
net$sel = ifelse((net$X == net$variable)==TRUE, paste("sel"), paste(""))
#or
net$sel = ifelse(net$X == net$variable, paste("sel"), paste(""))
net = subset(net, sel == "sel")
write.csv(net, file = "results/net_pseudo/net_sel.csv")
