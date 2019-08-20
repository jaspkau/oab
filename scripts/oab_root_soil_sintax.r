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
library(dendextend)
library(ggdendro)

# Make phyloseq object ----------------------------------------------------
#meta data sheet
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

d.s = merge_phyloseq(d.pc, d.pp)
d.ra = d.s

####Combine phyloseq objects

d.fin = merge_phyloseq(d_r, d.pc, d.pp)
d.fin

####relative abudances of phyla in roots and soil (Fig. X)

source("scripts/relative_abundance_plots.R")

###Filter d.fin for most abundant bacterial phyla observed in roots

d.fin2 = subset_taxa(d.fin, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")
d.fin2

####Sequence counts per sample
seq_counts = data.frame(sample_sums(d.fin2))
seq_counts = merge(seq_counts, sample_data(d), by = "row.names")
xlsx::write.xlsx(seq_counts, "results/seq_counts_pe.xlsx")

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d.fin2, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]
temp = subset(temp, Species == "P. praeclara")

# Once again, effective numbers to the resuce!
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef = 1/(1-temp$Simpson)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha = exp(temp$Shannon)

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$ef.sha)

alpha.kw = c()
for(i in c(14)){
  column = names(temp[i])
  k.demo = kruskal.test(ef.sha ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")
alpha.kw

avg = temp %>%
  group_by(Source) %>%
  summarise(simp = mean(ef.sha))
avg

# Hierarchial clustering --------------------------------------------------

sample_data(d.fin2)$int = paste(sample_data(d.fin2)$Species,".",
                                sample_data(d.fin2)$Source,".",
                                sample_data(d.fin2)$Pop_size,".",
                                sample_data(d.fin2)$replicate)
sample_data(d.fin2)$int = gsub(" ", "", sample_data(d.fin2)$int)
d2 = merge_samples(d.fin2, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int = vegdist(otu3, method = "bray")

otu3 = otu_table(as.matrix(otu3), taxa_are_rows = F)
d4 = merge_phyloseq(tax2, otu_table(as.matrix(otu3), 
                                    taxa_are_rows = F), sample_data(d2))
d.hc = tax_glom(d4, "Family")
tax.hc <- as(tax_table(d.hc),"matrix")
tax.hc <- as.data.frame(tax.hc)
otu.hm = merge(t(as.data.frame(otu_table(d.hc))), tax.hc, by = "row.names")
#otu.hm$rank = paste(as.character(otu.hm$Row.names),"|",substr(otu.hm$Family, 3, 5))
otu.hm$rank = paste(otu.hm$Family)
otu.hm = otu.hm[!(otu.hm$rank == "unidentified"),]
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))

who = names(sort(colMeans(otu.hm), decreasing = TRUE))[1:30]
otu.hm2 = otu.hm[,names(otu.hm) %in% who]

#weighted distance analysis
h = hclust(dist_w_int, method = "average")

dhc <- as.dendrogram(h) %>% set("labels_cex", 0.5)
ggd1 <- as.ggdend(dhc)

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g3 = heatmap.2(as.matrix(otu.hm2), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100),
               trace = "none")