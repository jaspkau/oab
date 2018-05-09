# Realtive abundance plots at OTU level ------------------------------------------------

d_f = merge_samples(d.fin2, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
gen_f$rank = paste(as.character(gen_f$Row.names),"|",substr(gen_f$Family, 3, 5))
list = as.character(gen_f$rank)
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank", "otu", "pval", "p.ad")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
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

p.otus = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.8, keyheight = 0.8))+ scale_y_continuous(labels = percent_format())
p.otus$data$variable = factor(p.otus$data$variable, ordered = TRUE, levels = rev(who))
p.otus

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d.fin2, taxrank = "Family")
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
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 3, name = "Pastel1"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 8, name = "Accent")))

library(scales)

p.fam = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.8, keyheight = 0.8))+ scale_y_continuous(labels = percent_format())
p.fam$data$variable = factor(p.fam$data$variable, ordered = TRUE, levels = rev(who))
p

#ggsave(file="jc.treatment.nms.jpg")

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d.fin2, taxrank = "Phylum")
d_f = merge_samples(d_f, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
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
who = names(sort(colMeans(gen_f), decreasing = TRUE))
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("blue2", "black","yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 8, name = "Accent")))

library(scales)

p.phy = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 1, reverse=T, keywidth = 0.8, keyheight = 0.8))+ scale_y_continuous(labels = percent_format())
p.phy$data$variable = factor(p.phy$data$variable, ordered = TRUE, levels = rev(who))



