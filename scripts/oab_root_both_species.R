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

###ROOT OAB ANALYSIS......................................
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "root")
decon
#decon = subset_samples(decon, Population2 != "X")
#decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence
source("scripts/decontaminate_phyloseq.R")

d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r
###Filter d.fin for most abundant bacterial phyla observed in roots

d_r = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")

d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

####OTU shared between species and unique ones

d.pc = subset_samples(d_r, Species == "P. cooperi")
d.pc
d.pc = prune_taxa(taxa_sums(d.pc) >= 1, d.pc)
d.pc
d.pp = subset_samples(d_r, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp
d.shared = prune_taxa(taxa_names(d.pp) %in% taxa_names(d.pc), d.pp)
d.shared

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d_r, measures = c("Shannon", "InvSimpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp <- ggplot(temp, aes(x=Species, y=InvSimpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (D)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$Simpson)

alpha.kw = c()
for(i in c(5)){
  column = names(temp[i])
  k.demo = kruskal.test(InvSimpson ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")

library(dplyr)

avg = temp %>%
  group_by(Species) %>%
  summarise(simpson = mean(InvSimpson))
avg

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = data.frame(otu_table(d_r))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

##PERMANOVA

###Weighted distance

a = adonis2(dist_w ~ sample_data(d3)$Species + sample_data(d3)$Pop_size + sample_data(d3)$Stage, 
            strata = sample_data(d3)$Species, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

# PCoA with bray --------------------------------------------------------------------

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
  geom_point(aes(color = Species), size=2) + state_col_ord +
  labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
       y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))
p

# Hierarchial clustering --------------------------------------------------

sample_data(d_r)$int = paste(sample_data(d_r)$Population,".",gsub("20","",sample_data(d_r)$Year))
sample_data(d_r)$int = gsub(" ", "", sample_data(d_r)$int)
d2 = merge_samples(d_r, "int")
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

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(5, 5), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

d.an = d.pc
d.an = tax_glom(d.pc, taxrank = "Family")

otu.hm = merge(t(as.data.frame(otu_table(d.an))), tax_table(d.an), by = "row.names")
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
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Pop_size", names(ancom.otu))]

anc = ANCOM(ancom.fin, multcorr = 2, sig = 0.05)
anc$detected
plot_ancom(anc)

library(devtools)  # Load the devtools package
#biocLite("DESeq2")
library(DESeq2)
#biocLite("WGCNA")
library(WGCNA)
#install_github("umerijaz/microbiomeSeq")  # Install the package
library(microbiomeSeq)  #load the package
#biocLite("GO.db")
library(GO.db)
#biocLite("impute")
library(impute)
#biocLite("preprocessCore")
library(preprocessCore)
source("http://www.bioconductor.org/biocLite.R")

####DeSeq
d_fam = taxa_level(d.pp, "Family")
#d_r = t(d_r)
##in the following command, make sure grouping label is stored as factor in met file
deseq_sig <- differential_abundance(d_fam, grouping_column = "Land_mg", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
deseq_sig$importance
p <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p)