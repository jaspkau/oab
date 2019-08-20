#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
#setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/oab/oab")

#https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/library("ape")
library("Biostrings")
library("ggplot2")
#devtools::install_github("GuangchuangYu/treeio")
#23

library(treeio)
library(ape)
library(ggtree)

#https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
#https://guangchuangyu.github.io/presentation/2016-ggtree-chinar/

###import associated matrix

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

#########read phyloseq object
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

taxa_names(d.fin) = gsub("otu", "denovo", taxa_names(d.fin))

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source)
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Stage)
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Year)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

##merge samples
dpop = merge_samples(d.fin, "int")

####for abudance based associated matrix
otu = data.frame(t(otu_table(dpop)))
row.names(otu) = gsub("otu", "denovo", row.names(otu))

otu2 = as.data.frame(t(otu))
otu2 = as.data.frame(t(decostand(otu2, method = "hellinger")))

###for raxml 

x = read.raxml("results/raxml_bacillales/RAxML_bipartitionsBranchLabels.bootFinal")

tree = ggtree(x, color="black", size=0.7) + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = bootstrap), vjust=-.5, hjust = 0.5, size=3)

tree$data$bootstrap = ifelse(tree$data$bootstrap < 50, paste(""), tree$data$bootstrap)

g = gheatmap(tree, otu2, offset = 0.0, width=0.5, font.size=3, colnames_angle=-45, hjust = 0, 
             low = "grey97", high = "black")
g

##make binary matrix
otu2 = as.data.frame(ifelse(otu == 0, 0, 1))
###

g = gheatmap(tree, otu2, offset = 0.0, width=0.5, font.size=3, colnames_angle=-45, hjust = 0, 
             low = "red3", high = "limegreen")
g

###for mrbayes
x = read.mrbayes("results/raxml_bacillales//result.nexusi.con.tre")
#x = read.raxml("results/phylo/tul_raxml/RAxML_bipartitionsBranchLabels.bootFinal")

tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=prob, vjust=-.5, size=3))

tree$data$prob[is.na(tree$data$prob)] = ""

tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, size=1))


##missleneaous

ggtree(x) +
  geom_label(mapping = aes(label = node), size = 2)

exp = ggtree(x, aes(color=branch.length), size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, hjust = 0.8, size=0.5)) +
  theme(legend.position="bottom")

cp = collapse(exp, node = 183)
cp + geom_point2(aes(subset=(node == 183)), size=5, shape=23, fill="steelblue")

####Collapse nodes

#Branch Length = consider an alignment of length 100, and RAxML 
#estimates a particular branch length to be 0.1,  
#the total number of substitutions on that branch should be 0.1 X 100
# collapse if prob < 0.6 and branch legnth is < 0.00

nod.claps = ifelse((tree$data$branch.length < 0.03)==TRUE, paste(tree$data$node), paste(""))
exp %>% collapse(nod.claps)

