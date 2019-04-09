aldiv = estimate_richness(d.fin2, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp1 <- ggplot(temp, aes(x=Source, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp1

aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp1 <- ggplot(temp, aes(x=Species, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp1

d.r.pc = subset_samples(d_r, Species == "P. cooperi")
d.r.pc = prune_taxa(taxa_sums(d.r.pc) >= 1, d.r.pc)

aldiv = estimate_richness(d.r.pc, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp1 <- ggplot(temp, aes(x=Pop_size, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp1

bp2 <- ggplot(temp, aes(x=Stage, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp2

bp3 <- ggplot(temp, aes(x=Year, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp3


d.r.pc = subset_samples(d_r, Species == "P. praeclara")
d.r.pc = prune_taxa(taxa_sums(d.r.pc) >= 1, d.r.pc)

aldiv = estimate_richness(d.r.pc, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

bp4 <- ggplot(temp, aes(x=Pop_size, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp4

bp5 <- ggplot(temp, aes(x=Stage, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp5

bp6 <- ggplot(temp, aes(x=Year, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp6

library(ggpubr)
ggarrange(bp1, bp2, bp3, bp4, bp5, bp6, nrow = 2, ncol = 3)
