####decontaminate phyloseq object based on frequency and prevelence
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

df <- as.data.frame(sample_data(decon)) # Put sample_data into a ggplot-friendly d
df$LibrarySize <- sample_sums(decon)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
ggsave(p, width = 8, height = 6, units = "in", file="results/lib_size.jpg")

###combined method of decomtamination based on prevelanec and frequency
sample_data(decon)$is.neg <- sample_data(decon)$Sample_or_Control == "Control"
table(get_variable(decon, "is.neg"), useNA="always")
#contamdf.prev <- isContaminant(decon, method="combined", neg="is.neg", conc="DNA_conc", threshold=0.1)
#contamdf.prev <- isContaminant(decon, method="prevalence", neg="is.neg", threshold=0.1)
contamdf.prev <- isContaminant(decon, method="frequency", neg="is.neg", conc="DNA_conc", threshold=0.1)
#contamdf.prev <- isContaminant(decon, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon.d <- prune_taxa(!contamdf.prev$contaminant, decon)
decon.d