otu <- read.delim(file = "data/otu_table_no_singletons_sintax.txt", 
                  sep = "\t", header = T)
otu = otu[,-ncol(otu)]
row.names(otu) = paste(gsub("denovo", "o", otu[,1]))
otu = otu[,-1]
#Rarefy(otu, depth = min(rowSums(otu)))
otu = otu[,colSums(otu) > 0]
site_list = colnames(otu)
otu = t(otu)
otu_tab = otu_table(as.matrix(otu), taxa_are_rows = F)

###Format SINTAX taxonomy table

library(reshape2)

tax = read.delim(file = "data/tax.sintax", sep = "\t", header = F)
row.names(tax) = tax$V1
list = tax$V2
tax2 = colsplit(list, pattern ="\\(|\\),", c("Kingdom", "Kingdom_conf", "Phylum", "Phylum_conf", "Class", "Class_conf", "Order", "Order_conf", "Family", "Family_conf", "Genus", "Genus_conf", "Species", "Species_conf"))
source("scripts/correct_tax_dataframe.R")
tax2$Species_conf = gsub("\\)", "", tax2$Species_conf)
tax2$Species_conf = as.numeric(tax2$Species_conf)
tax2[is.na(tax2)] <- 0
row.names(tax2) = row.names(tax)

source("scripts/tax_func.R") #90% conf
level = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
for (i in 1:7){ ###i 1 = kingdom, i 2 = phylum etc
  tax2 = tax_assign(tax2, paste(level[i], "_conf", sep = ""), level[i])
}

tax2 = tax2[,-c(2,4,6,8,10,12,14)]

tax2$row = row.names(tax2)
tax2[,8:9] = colsplit(tax2$row, " ", c("otu", "seq"))
row.names(tax2) = paste(gsub("denovo", "o", tax2$row))
tax2 = tax2[,-c(8:9)]
tax2 = tax_table(as.matrix(tax2))

#meta data
row.names(met) = met$code
met$Population = as.factor(met$Population)
met$Species = as.factor(met$Species)
met$Year = as.factor(met$Year)
met$Stage = as.factor(met$Stage)
met$Samp_L = as.factor(met$Samp_L)
met$Pop_size = as.factor(met$Pop_size)
met$Sample_or_Control = as.factor(met$Sample_or_Control)
met$pop.year = paste(met$Population, ".", gsub("20", "", met$Year))
met$pop.year = gsub(" ", "", met$pop.year)
met$pop.sta = paste(met$Population, ".", met$Stage)
met$pop.sta = gsub(" ", "", met$pop.sta)
met$int = paste(met$Population, ".", met$Stage, ".", gsub("20", "", met$Year))
met$int = gsub(" ", "", met$int)

#phyloseq object

d = merge_phyloseq(tax2, otu_tab, sample_data(met))
d
d = subset_taxa(d, Kingdom == "d:Bacteria")
d
d = subset_taxa(d, Phylum != "p:Cyanobacteria/Chloroplast")
d  
