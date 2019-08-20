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

###ROOT OMF ANALYSIS......................................
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

#phyloseq object

decon = subset_samples(d, Source == "root")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

decon.d = subset_samples(decon.d, Species == "P. praeclara")
decon.d
#decon.d = subset_samples(decon.d, Population2 != "X")
#decon.d
d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

###Filter for most abundant bacterial phyla observed in roots

d_r = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")

d_r

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef = 1/(1-temp$Simpson)
# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha = exp(temp$Shannon)

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$ef.sha)

alpha.kw = c()
for(i in c(12)){
  column = names(temp[i])
  k.demo = kruskal.test(ef.sha ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")
alpha.kw

avg = temp %>%
  group_by(Pop_size) %>%
  summarise(simp = mean(ef.sha))
avg

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

d.an = d_r
d.an = tax_glom(d_r, taxrank = "Family")
tax.an <- as(tax_table(d.an),"matrix")
tax.an <- as.data.frame(tax.an)
otu.hm = merge(t(as.data.frame(otu_table(d.an))), tax.an, by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
otu.hm$rank = paste(as.character(otu.hm$Row.names),"|",substr(otu.hm$Family, 3, 5))
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))  ##columns = OTUs and should be counts

ancom.otu = merge(otu.hm, sample_data(d.an), by = "row.names")
row.names(ancom.otu) = ancom.otu$Code
ancom.otu = ancom.otu[,-1]
names(ancom.otu)
##look for the grouping variable you want to use
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Stage", names(ancom.otu))]

library(ancom.R)
anc = ANCOM(ancom.fin, sig = 0.05)
anc$detected
plot_ancom(anc)