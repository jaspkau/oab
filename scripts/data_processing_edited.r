library(phyloseq)
library(ggplot2)
library(reshape2)
library(permute)
library(lattice)
library(vegan)
library(plyr)

setwd("C:/Users/js/Desktop/R24april/New folder/")

######
#Format the taxonomy table
#####
tax = read.delim(file = "tax.txt", sep = "\t", header = F)
row.names(tax) = tax$V1
list = tax$V2
tax2 = colsplit(list, pattern ="\\(|\\),", c("Kingdom", "Kingdom_conf", "Phylum", "Phylum_conf", "Class", "Class_conf", "Order", "Order_conf", "Family", "Family_conf", "Genus", "Genus_conf"))
tax2$Genus_conf = gsub("\\)", "", tax2$Genus_conf)
row.names(tax2) = row.names(tax)
tax2$Genus_conf = as.numeric(tax2$Genus_conf)


tax_assign = function(x, conf, rank){
  list = x[conf][[conf]]
  list2 = as.character(x[rank][[rank]])
  list3 = ifelse(list < 95, "unknown", list2)
  x[rank] = list3
  return(x)
}


level = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
for (i in 1:6){
  tax2 = tax_assign(tax2, paste(level[i], "_conf", sep = ""), level[i])
}

tax2 = tax2[,-c(2,4,6,8,10,12)]
tax2 = tax_table(as.matrix(tax2))  #make the taxonomy table as a phyloseq

#####
#Format the otu table
#####
otu = read.delim(file = "otu_table.txt", sep = "\t", header = T)
row.names(otu) = otu[,1]
otu = otu[,-1]
names(otu) = gsub("_", ".", names(otu))
otu = otu_table(as.matrix(otu), taxa_are_rows = T)

#####
#Bring in the tree
#####
#http://www.phytools.org/read.newick/v0.5/read.newick.R
tre = read.newick(file = "tree.tre")
tre$tip.label = gsub("_", ":", tre$tip.label)
tre$tip.label = gsub(":MS28F", "_MS28F", tre$tip.label)

#####
#Set up the metadata
#####
met = read.csv(file = "met.csv", header = T)
row.names(met) = paste(met$Sample, ".MS28F", sep = "")

#####
#Create the phyloseq object
#####
d = merge_phyloseq(tax2, otu, tre, sample_data(met))


d <- prune_taxa(taxa_sums(d) > 0, d)


#####
#Calculate unifrac
#####
unif_w = UniFrac(d, weighted=T, normalized = T)
unif_uw = UniFrac(d, weighted=F, parallel=T)

unif_w = as.matrix(unif_w)
unif_uw = as.matrix(unif_uw)

wuf = as.dist(unif_w)
uwuf = as.dist(unif_uw)

#####
#alpha diversity
#####
r = estimate_richness(d)
r = merge(r, met, by = "row.names")

plot(r$Group, r$Chao1)
#####
#beta diversity
#####
a = adonis(wuf ~ sample_data(d)$Group)

pc = capscale(uwuf~1)
cap = data.frame(pc$CA$u)
cap = merge(met, cap, by = "row.names")
rownames(cap) = cap[,1]
cap = cap[,2:ncol(cap)]
s = summary(pc)
p = ggplot(cap, aes(MDS1, MDS2)) + geom_point(aes(color = Group), size = 5) + theme_bw(base_size = 20)
p = p + labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
             y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(p$data, c("Group"), find_hull)
p = p + geom_polygon(data = hulls, aes(fill = Group), alpha = 0.15, linetype = 0)



