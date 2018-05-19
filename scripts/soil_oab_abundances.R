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

###ROOT OMF ANALYSIS......................................
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

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

r.otus = taxa_names(d_r)

###make soil phyloseq

decon = subset_samples(d, Source == "soil")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence
<<<<<<< HEAD

sample_data(decon)$is.neg <- sample_data(decon)$Sample_or_Control == "Control"
#contamdf.prev <- isContaminant(decon, method="frequency", neg="is.neg", conc="DNA_conc", threshold=0.5)
contamdf.prev <- isContaminant(decon, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon.d <- prune_taxa(!contamdf.prev$contaminant, decon)
decon.d
=======
source("scripts/decontaminate_phyloseq.R")
>>>>>>> 6a1b8b269f85ee301136e795be245bcf9c854969

d_s = subset_samples(decon.d, Sample_or_Control == "Sample")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
<<<<<<< HEAD
d_s = merge_samples(d_s, "Population")

d_oab = prune_taxa(r.otus, d_s) ##extract root otus from soil phyloseq
alltaxa = taxa_names(d_s)
nonoab <- alltaxa[!(alltaxa %in% r.otus)]
d_nonoab = prune_taxa(nonoab, d_s)

oab.sum = sample_sums(d_oab)
nonoab.sum = sample_sums(d_nonoab)

abun = cbind(oab.sum, nonoab.sum)

# Realtive abundance plots at OTU level ------------------------------------------------

gen_f = data.frame(abun/rowSums(abun))
who = names(gen_f)
gen_f$sl = row.names(gen_f)
m = melt(gen_f, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("gray", "black"))
                                                         
library(scales)

p.otus = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.8, keyheight = 0.8))+ scale_y_continuous(labels = percent_format())
p.otus$data$variable = factor(p.otus$data$variable, ordered = TRUE, levels = rev(who))
p.otus
=======
d_s

otu_s = data.frame(otu_table(d_s))
otu_s2 = otu_s[r.otus,]
otu_s2 = otu_s2[,colSums(otu_s2) > 0]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_oab = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_oab
d_oab = prune_taxa(taxa_sums(d_oab) >= 1, d_s)
d_oab
>>>>>>> 6a1b8b269f85ee301136e795be245bcf9c854969
