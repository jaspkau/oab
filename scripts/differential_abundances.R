#source("https://bioconductor.org/biocLite.R")
library(devtools)  # Load the devtools package
#biocLite("DESeq2")
library(DESeq2)
#install_github("umerijaz/microbiomeSeq")  # Install the package
library(microbiomeSeq)  #load the package
#biocLite("GO.db")
library(GO.db)
#biocLite("impute")
library(impute)
#biocLite("preprocessCore")
library(preprocessCore)
source("http://www.bioconductor.org/biocLite.R")

d_sps = merge_samples(d_r, "Species")
#install.packages("RAM")
library(RAM)

share_df = data.frame(otu_table(d_sps))
share_df$taxonomy = seq(1:nrow(share_df))
shared.otus = shared.OTU(share_df)
shared <- shared.OTU(data=list(share_df=share_df))

####DeSeq
d_fam = taxa_level(d_r, "Family")
#d_r = t(d_r)
##in the following command, make sure grouping label is stored as factor in met file
deseq_sig <- differential_abundance(d_fam, grouping_column = "Stage", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
p <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p)

####Co-occurance
physeq <- taxa_level(d_r, "Phylum")
co_occr <- co_occurence_network(physeq, grouping_column = "Pop_size", rhos = 0.20, select.condition = "S", scale.vertex.size=3, scale.edge.width=15)

