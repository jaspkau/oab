setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
setwd("C:/Users/jas/Google Drive/Metagenomics/oab/")
setwd("/Users/administrator/Desktop/oab/")

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

###soil OMF ANALYSIS......................................
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
decon.d = subset_samples(decon.d, Population2 != "X")
decon.d
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
d.pp = subset_samples(d.pp, Population2 != "X")
d.pp
d.pp = subset_samples(d.pp, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp

####Combine phyloseq objects

d.fin = merge_phyloseq(d_r, d.pp)
d.fin

###Filter d.fin for most abundant bacterial phyla observed in roots

d.fin2 = subset_taxa(d.fin, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria"| Phylum == "p:Bacteroidetes")
d.fin2

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d.fin2, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp <- ggplot(temp, aes(x=Source, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$Simpson)

alpha.kw = c()
for(i in c(14)){
  column = names(temp[i])
  k.demo = kruskal.test(Simpson ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = data.frame(otu_table(d.fin2))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA
###Weighted distance

a = adonis2(dist_w ~ sample_data(d3)$Source, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different

sample_data(d.fin2)$int = paste(sample_data(d.fin2)$Source,".",sample_data(d.fin2)$Population,".",gsub("20","",sample_data(d.fin2)$Year))
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

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

ancom.otu = data.frame(otu_table(d.fin2)) ##columns = OTUs and should be counts
ancom.otu = merge(ancom.otu, sample_data(d.an), by = "row.names")
row.names(ancom.otu) = ancom.otu$Code
ancom.otu = ancom.otu[,-1]
names(ancom.otu)
##look for the grouping variable you want to use
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Source", names(ancom.otu))]

anc = ANCOM(ancom.fin, multcorr = 1, sig = 0.05)
anc$detected
plot_ancom(anc)

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(5, 5), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")