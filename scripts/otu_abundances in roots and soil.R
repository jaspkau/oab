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

#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Species == "P. cooperi")
decon
decon = subset_samples(decon, Sample_or_Control == "Sample")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon
d_r = subset_samples(decon, Source == "root")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

d.abun = merge_samples(decon, "Source")
mytaxa = taxa_names(d_r)
d.abun2 = prune_taxa(mytaxa, d.abun)

otu3 = data.frame(otu_table(d.abun2))
otu3 = decostand(otu3, method = "hellinger")
otu3 = otu_table(otu3, taxa_are_rows = FALSE)
d.abun3 = merge_phyloseq(otu3, tax_table(d.abun2), sample_data(d.abun2))
test = data.frame(otu_table(d.abun3))
test2 = data.frame(t(test))
plot(test2)
label = row.names(test2)

p = ggplot(test2, aes(root, soil)) + geom_point() + 
  geom_text(aes(label=label),size = 3, vjust = "inward", 
            hjust = "inward", check_overlap = TRUE)
p

# PCoA with bray --------------------------------------------------------------------

d.pcoa = prune_taxa(mytaxa, decon)
library(vegan)
otu2 = data.frame(otu_table(d.pcoa))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

state_col_ord = scale_color_manual(values=c("red", "black"))

####Weighted

pc = capscale(dist_w ~ 1, comm = rel_otu_code) ###~ means function of nothing
pc$CA$eig
s = summary(pc)
cap = data.frame(pc$CA$u)
plot(cap[,1], cap[,2])
cap = merge(sample_data(d), cap, by = "row.names")
#cap$Population = cap$Row.names
label = cap$int

p = ggplot(cap, aes(x= MDS1, y= MDS2))+theme_bw(base_size = 15) +
  geom_point(aes(color = Source), size=2) + state_col_ord +
  labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
       y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))
p

                  