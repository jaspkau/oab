
# Realtive abundance plots at Phylum level ------------------------------------------------

d_f = tax_glom(d.fin, taxrank = "Phylum")
sample_data(d_f)$int = paste(sample_data(d_f)$Species,".", sample_data(d_f)$Source)
d_f = merge_samples(d_f, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
tax.df <- as(tax_table(d_f),"matrix")
tax.df <- as.data.frame(tax.df)
gen_f = merge(gen_f, tax.df, by = "row.names")
gen_f$rank = as.character(gen_f$Phylum)
#gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:10]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)

state_col2 = scale_fill_manual(values = c(brewer.pal(n = 8, name = "Dark2"),
                                                brewer.pal(n = 11, name = "Spectral"), "blue2", "black","yellow1",
                                                "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                "slategray4", "seagreen4" , "aquamarine",
                                                "tomato2"))

library(scales)

p.phy = ggplot(m, aes(sl, fill = variable)) + 
  geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + 
  state_col2 + xlab("Sample") + 
  ylab("Relative Abundance (%)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(axis.text.y = element_text(hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "plain", size = 12, family = "Arial")) + 
  guides(fill = guide_legend(ncol = 1, reverse=FALSE, keywidth = 0.8, keyheight = 0.8)) + 
  scale_y_continuous(labels = percent_format())
p.phy$data$variable = factor(p.phy$data$variable, ordered = TRUE, levels = who)

p.phy

#ggsave(file="jc.treatment.nms.jpg")