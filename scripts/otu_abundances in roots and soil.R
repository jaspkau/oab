setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
setwd("C:/Users/jaspr/Google Drive/Metagenomics/oab/")
setwd("/Users/administrator/Documents/jaspreet/oab/oab")

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "root")
decon

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

d.fin2 = subset_taxa(d.fin, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria"| Phylum == "p:Bacteroidetes")

d.fin2 = subset_taxa(d.fin, Family == "f:Pseudomonadaceae")

d.abun = merge_samples(d.fin2, "Source")
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

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

ancom.otu = t(data.frame(otu_table(d.an))) ##columns = OTUs and should be counts
ancom.otu = merge(ancom.otu, sample_data(d.an), by = "row.names")
row.names(ancom.otu) = ancom.otu$Code
ancom.otu = ancom.otu[,-1]
names(ancom.otu)
##look for the grouping variable you want to use
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Pop_size", names(ancom.otu))]

anc = ANCOM(ancom.fin, multcorr = 1, sig = 0.05)
anc$detected
plot_ancom(anc)
