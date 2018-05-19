
# Correct the taxonomy database -------------------------------------------

###fix order assignments

tax2$Order = ifelse(grepl('o:',tax2$Class)==TRUE, tax2$Class, tax2$Order)
tax2$Order_conf = ifelse(grepl('o:',tax2$Class)==TRUE, tax2$Class_conf, tax2$Order_conf)
tax2$Class_conf = ifelse(grepl('o:',tax2$Class)==TRUE, paste(""), tax2$Class_conf)
tax2$Class = ifelse(grepl('o:',tax2$Class)==TRUE, paste("unidentified"), tax2$Class)

tax2$Family = ifelse(grepl('f:',tax2$Class)==TRUE, tax2$Class, tax2$Family)
tax2$Family_conf = ifelse(grepl('f:',tax2$Class)==TRUE, tax2$Class_conf, tax2$Family_conf)
tax2$Class_conf = ifelse(grepl('f:',tax2$Class)==TRUE, paste(""), tax2$Class_conf)
tax2$Class = ifelse(grepl('f:',tax2$Class)==TRUE, paste("unidentified"), tax2$Class)

tax2$Genus = ifelse(grepl('g:',tax2$Class)==TRUE, tax2$Class, tax2$Genus)
tax2$Genus_conf = ifelse(grepl('g:',tax2$Class)==TRUE, tax2$Class_conf, tax2$Genus_conf)
tax2$Class_conf = ifelse(grepl('g:',tax2$Class)==TRUE, paste(""), tax2$Class_conf)
tax2$Class = ifelse(grepl('g:',tax2$Class)==TRUE, paste("unidentified"), tax2$Class)

tax2$Species = ifelse(grepl('s:',tax2$Class)==TRUE, tax2$Class, tax2$Species)
tax2$Species_conf = ifelse(grepl('s:',tax2$Class)==TRUE, tax2$Class_conf, tax2$Species_conf)
tax2$Class_conf = ifelse(grepl('s:',tax2$Class)==TRUE, paste(""), tax2$Class_conf)
tax2$Class = ifelse(grepl('s:',tax2$Class)==TRUE, paste("unidentified"), tax2$Class)

####fix Order assignment

tax2$Family = ifelse(grepl('f:',tax2$Order)==TRUE, tax2$Order, tax2$Family)
tax2$Family_conf = ifelse(grepl('f:',tax2$Order)==TRUE, tax2$Order_conf, tax2$Family_conf)
tax2$Order_conf = ifelse(grepl('f:',tax2$Order)==TRUE, paste(""), tax2$Order_conf)
tax2$Order = ifelse(grepl('f:',tax2$Order)==TRUE, paste("unidentified"), tax2$Order)

tax2$Genus = ifelse(grepl('g:',tax2$Order)==TRUE, tax2$Order, tax2$Genus)
tax2$Genus_conf = ifelse(grepl('g:',tax2$Order)==TRUE, tax2$Order_conf, tax2$Genus_conf)
tax2$Order_conf = ifelse(grepl('g:',tax2$Order)==TRUE, paste(""), tax2$Order_conf)
tax2$Order = ifelse(grepl('g:',tax2$Order)==TRUE, paste("unidentified"), tax2$Order)

tax2$Species = ifelse(grepl('s:',tax2$Order)==TRUE, tax2$Order, tax2$Species)
tax2$Species_conf = ifelse(grepl('s:',tax2$Order)==TRUE, tax2$Order_conf, tax2$Species_conf)
tax2$Order_conf = ifelse(grepl('s:',tax2$Order)==TRUE, paste(""), tax2$Order_conf)
tax2$Order = ifelse(grepl('s:',tax2$Order)==TRUE, paste("unidentified"), tax2$Order)

#####fix family assigment

tax2$Genus = ifelse(grepl('g:',tax2$Family)==TRUE, tax2$Family, tax2$Genus)
tax2$Genus_conf = ifelse(grepl('g:',tax2$Family)==TRUE, tax2$Family_conf, tax2$Genus_conf)
tax2$Family_conf = ifelse(grepl('g:',tax2$Family)==TRUE, paste(""), tax2$Family_conf)
tax2$Family = ifelse(grepl('g:',tax2$Family)==TRUE, paste("unidentified"), tax2$Family)

tax2$Species = ifelse(grepl('s:',tax2$Family)==TRUE, tax2$Family, tax2$Species)
tax2$Species_conf = ifelse(grepl('s:',tax2$Family)==TRUE, tax2$Family_conf, tax2$Species_conf)
tax2$Family_conf = ifelse(grepl('s:',tax2$Family)==TRUE, paste(""), tax2$Family_conf)
tax2$Family = ifelse(grepl('s:',tax2$Family)==TRUE, paste("unidentified"), tax2$Family)

####fix genus assignment

tax2$Species = ifelse(grepl('s:',tax2$Genus)==TRUE, tax2$Genus, tax2$Species)
tax2$Species_conf = ifelse(grepl('s:',tax2$Genus)==TRUE, tax2$Genus_conf, tax2$Species_conf)
tax2$Genus_conf = ifelse(grepl('s:',tax2$Genus)==TRUE, paste(""), tax2$Genus_conf)
tax2$Genus = ifelse(grepl('s:',tax2$Genus)==TRUE, paste("unidentified"), tax2$Genus)

