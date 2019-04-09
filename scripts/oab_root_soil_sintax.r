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

###Filter d.fin for most abundant bacterial phyla observed in roots

d.fin2 = subset_taxa(d.fin, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                       Phylum == "p:Actinobacteria")
d.fin2

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d.s, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

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

# Beta diversity with bray ------------------------------------------------
library(vegan)

####model 1, root com ~ species, and population size and stage

d.mod1 = subset_samples(d.fin2, Source == "root")
d.mod1
d.mod1 = prune_taxa(taxa_sums(d.mod1) >= 1, d.mod1)
d.mod1

otu2 = data.frame(otu_table(d_r))
otu2 = decostand(otu2, method = "hellinger")
#rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d.fin2))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

##PERMANOVA

###Weighted distance
a = adonis2(dist_w ~ sample_data(d3)$Species +
              sample_data(d3)$Pop_size +
              sample_data(d3)$Stage, 
            starata = sample_data(d3)$Species, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

############################model 2, soil com ~ species + popsize

d.mod2 = subset_samples(d.fin2, Source == "soil")
d.mod2
d.mod2 = prune_taxa(taxa_sums(d.mod2) >= 1, d.mod2)
d.mod2

otu2 = data.frame(otu_table(d.mod2))
otu2 = decostand(otu2, method = "hellinger")
#rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d.mod2))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

a = adonis2(dist_w ~ sample_data(d3)$Species +
              sample_data(d3)$Pop_size, 
            strata = sample_data(d3)$Species, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

########################model 3 = root and soil diversity in four vegetation management

d.mod3 = subset_samples(d.fin2, Population == "MNP")
d.mod3
d.mod3 = prune_taxa(taxa_sums(d.mod3) >= 1, d.mod3)
d.mod3

otu2 = data.frame(otu_table(d.mod3))
otu2 = decostand(otu2, method = "hellinger")
#rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d.mod3))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

a = adonis2(dist_w ~ sample_data(d3)$Land_mg, strata = sample_data(d3)$Source, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

####model 4 = soil diversity between orchid occupied and unoccupied locations

d.mod4 = subset_samples(d.fin2, Species == "P. praeclara")
d.mod4
d.mod4 = subset_samples(d.mod4, Source == "soil")
d.mod4
d.mod4 = prune_taxa(taxa_sums(d.mod4) >= 1, d.mod4)
d.mod4

otu2 = data.frame(otu_table(d.mod4))
otu2 = decostand(otu2, method = "hellinger")
#rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d.mod4))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

a = adonis2(dist_w ~ sample_data(d3)$Samp_L, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

# Hierarchial clustering --------------------------------------------------

sample_data(d.mod3)$int = paste(sample_data(d.mod3)$Source,".",sample_data(d.mod3)$Land_mg)
sample_data(d.mod3)$int = gsub(" ", "", sample_data(d.mod3)$int)
d2 = merge_samples(d.mod3, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int = vegdist(otu3, method = "bray")

otu3 = otu_table(as.matrix(otu3), taxa_are_rows = F)
d4 = merge_phyloseq(tax2, otu_table(as.matrix(otu3), 
                                    taxa_are_rows = F), sample_data(d2))
d.hc = tax_glom(d4, "Family")

otu.hm = merge(t(as.data.frame(otu_table(d.hc))), tax_table(d.hc), by = "row.names")
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

p1 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p1

p2 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p2

p3 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p3

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g3 = heatmap.2(as.matrix(otu.hm2), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100),
               trace = "none")

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

d.pp2 = subset_taxa(d.pp, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                               Phylum == "p:Actinobacteria")
d.pp2
d.an = d.pp2
d.an = tax_glom(d.pp2, taxrank = "Family")

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
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Samp_L", names(ancom.otu))]

anc = ANCOM(ancom.fin, multcorr = 2, sig = 0.05)
anc$detected
plot_ancom(anc)