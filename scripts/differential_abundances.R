setwd("/Users/administrator/Documents/jaspreet/oab/oab/")

#source("https://bioconductor.org/biocLite.R")
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
decon = subset_samples(decon, Population2 != "X")
decon

####decontaminate phyloseq object based on frequency and prevelence
source("scripts/decontaminate_phyloseq.R")

d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
###Filter d.fin for most abundant bacterial phyla observed in roots

d_r = subset_taxa(d_r, Phylum == "p:Proteobacteria"| Phylum == "p:Firmicutes"|
                    Phylum == "p:Actinobacteria")

d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

d_sps = merge_samples(d_r, "Species")
#install.packages("RAM")
library(RAM)

share_df = data.frame(otu_table(d_sps))
share_df$taxonomy = seq(1:nrow(share_df))
shared.otus = shared.OTU(share_df)
shared <- shared.OTU(data=list(share_df=share_df))

####DeSeq
d_fam = taxa_level(d_r, "Order")
#d_r = t(d_r)
##in the following command, make sure grouping label is stored as factor in met file
deseq_sig <- differential_abundance(d_fam, grouping_column = "Species", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
p <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p)