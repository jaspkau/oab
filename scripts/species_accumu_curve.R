setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab")
setwd("C:/Users/jas/Google Drive/Metagenomics/oab")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/oab/oab")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(readxl)

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

# Species accumulation curve -------------------------------------------------------

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source,".",sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source,".",as.factor(sample_data(d.fin)$Year))
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

sample_data(d_r)$int = paste(sample_data(d_r)$Species,".",sample_data(d_r)$Stage)
sample_data(d_r)$int = gsub(" ", "", sample_data(d_r)$int)

dsac = merge_samples(d.fin, "int")

otu_rc = data.frame(otu_table(dsac)) ###rows should be samples

###
library(iNEXT)
otu_rc = data.frame(t(otu_rc)) ####columns should be samples
m <- c(100, 1000, 2000, 10000, 40000, 800000)
out = list(
  iNEXT(otu_rc, q=1, datatype="abundance", size=m, nboot = 1))
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g1 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g1

library(parallel)
cl <- makeCluster(20)
out <- clusterApply( 
  cl,
  iNEXT(otu_rc, q=1, datatype="abundance", size=m, nboot = 5),
  function(f) f()
)
stopCluster(cl)
