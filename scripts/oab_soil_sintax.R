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

# Make phyloseq object ----------------------------------------------------
#meta data sheet
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "soil")
decon
decon = subset_samples(decon, Species == "P. praeclara")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency and prevelence

sample_data(decon)$is.neg <- sample_data(decon)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(decon, method="prevalence", neg="is.neg", conc="DNA_conc", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon.d <- prune_taxa(!contamdf.prev$contaminant, decon)
decon.d

d_s = subset_samples(decon.d, Sample_or_Control == "Sample")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

# Alpha diversity ---------------------------------------------------------

temp = estimate_richness(d_s)
temp = merge(met, temp, by = "row.names")
p =  ggplot(temp, aes(Population, Chao1))+ geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

a = summary(aov(Shannon ~ Population + Year, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad
a = summary(aov(Simpson ~ Population + Year, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad

plot_richness(d_s, x= "Population", measures=c("Shannon", "Simpson") )

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d_s))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

a = adonis2(dist_w ~ sample_data(d3)$Population + as.factor(sample_data(d3)$Year), permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different

d2 = merge_samples(d_s, "pop.year")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
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
               Rowv = as.dendrogram(h), margins = c(10, 3), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")

# Realtive abundance plots at OTU level ------------------------------------------------

d_f = merge_samples(d_s, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = as.character(gen_f$Row.names)
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:50]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 5, name = "Blues"),brewer.pal(n = 10, name = "Paired"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "gold1", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 11, name = "Spectral")))

library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.8, keyheight = 0.8))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p
ggsave(p, width = 8, height = 6, units = "in", file="results/rel_abun_otu.jpg")

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_s, taxrank = "Family")
d_f = merge_samples(d_f, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = as.character(gen_f$Family)
#gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
gen_f$rank = ifelse(gen_f$Phylum == "unidentified", paste(as.character(gen_f$Kingdom), as.character(gen_f$Phylum), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Phylum != "unidentified" &  gen_f$Class == "unidentified", paste(as.character(gen_f$Phylum), as.character(gen_f$Class), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Class != "unidentified" &  gen_f$Order == "unidentified", paste(as.character(gen_f$Class), as.character(gen_f$Order), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Order != "unidentified" &  gen_f$Family == "unidentified", paste(as.character(gen_f$Order), as.character(gen_f$Family), sep = ";"), gen_f$rank)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "Family", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 6, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 6)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.5, keyheight = 0.5))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))

p

ggsave(p, width = 8, height = 6, units = "in", file="results/rel_abun_fam.jpg")

