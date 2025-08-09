#files made from nmds plot files...
metavallTMT=metamergedunit
OTUallTMT=alltmtmergedunitOTUmdsdf 
Tax=OTUClass_df.mat
View(Tax)
View(OTUallTMT)
View(metavallTMT)

###### Exporting Cleaned porifera_df file #####
# this will save the file to your working directory
row.names(alltmtmergedunitOTUmdsdf) <- alltmtmergedunitOTUmdsdf$UNIT #make Units rowname for OTU table and metadata file.

#Making otu file with new sample names as row names
alltmtmergedunitOTUmdsdf_UNIT<-alltmtmergedunitOTUmdsdf
names(alltmtmergedunitOTUmdsdf_UNIT)
alltmtmergedunitOTUmdsdf_UNIT<-alltmtmergedunitOTUmdsdf_UNIT[,-(1:3)]
names(alltmtmergedunitOTUmdsdf_UNIT)


#making metadata file as well. 
alltmtmergedunitOTUmdsdf_metadata<-alltmtmergedunitOTUmdsdf
names(alltmtmergedunitOTUmdsdf_metadata)
alltmtmergedunitOTUmdsdf_metadata<-alltmtmergedunitOTUmdsdf_metadata[,-(4:96)]
names(alltmtmergedunitOTUmdsdf_metadata)

#getting that data into a phyloseq object so that we can plot it
```{r setup, include=FALSE}
#Preparing files for phyloseq from cumulative diversity file. All data from all time points combined...
#Make dataframes matrices 
#making data file a matrix
alltmtmergedunitOTUmdsdf.mat<-data.matrix(alltmtmergedunitOTUmdsdf, rownames.force = NA)
row.names(alltmtmergedunitOTUmdsdf.mat)
class(alltmtmergedunitOTUmdsdf.mat)
dim(alltmtmergedunitOTUmdsdf.mat)

#checking classification file a matrices
View(alltmtmergedunitOTUmdsdf.mat)
View(taxmat)

otumat = alltmtmergedunitOTUmdsdf.mat
taxmat = OTUClass_df.mat
head(taxmat)
head(otumat)

sampledata=alltmtmergedunitOTUmdsdf_metadata
sampledata=sample_data(alltmtmergedunitOTUmdsdf_metadata)
sampledata

OTU = otu_table(otumat, taxa_are_rows = FALSE)
TAX =tax_table(taxmat)


cumporphyseq = phyloseq(OTU, TAX, sampledata)

cumporphyseq

ntaxa (cumporphyseq)
nsamples (cumporphyseq)
rank_names (cumporphyseq)

#Calculate number of species per each individual major sponge group...##
cumporphyseq_dem_nOTU=subset_taxa(cumporphyseq, Class=="Demospongiae")
#Prune OTUs that are not present in any of the samples
cumporphyseq_dem_nOTU=prune_taxa(taxa_sums(cumporphyseq_dem_nOTU) > 0, cumporphyseq_dem_nOTU) #This actually gets rid of samples that have more than 0 so not exactly sure what prune_taxa is doing. Need to figure this out... Made a heatmap with and without
ntaxa(cumporphyseq_dem_nOTU) #58 species for Demospongiea

#calcarea
cumporphyseq_calc_nOTU=subset_taxa(cumporphyseq, Class=="Calcarea")
#Prune OTUs that are not present in any of the samples
umporphyseq_calc_nOTU=prune_taxa(taxa_sums(cumporphyseq_calc_nOTU) > 0, cumporphyseq_calc_nOTU)
ntaxa(umporphyseq_calc_nOTU) #19 species for Calcarea

#Homoslceromorpha
cumporphyseq_homo_nOTU=subset_taxa(cumporphyseq, Class=="Homoscleromorpha")
#Prune OTUs that are not present in any of the samples
umporphyseq_homo_nOTU=prune_taxa(taxa_sums(cumporphyseq_homo_nOTU) > 0, cumporphyseq_homo_nOTU)
ntaxa(umporphyseq_homo_nOTU) #8 species for Homoscleromorpha

#Keratosa
cumporphyseq_ker_nOTU=subset_taxa(cumporphyseq, Subclass=="Keratosa")
#Prune OTUs that are not present in any of the samples
umporphyseq_ker_nOTU=prune_taxa(taxa_sums(cumporphyseq_ker_nOTU) > 0, cumporphyseq_ker_nOTU)
ntaxa(umporphyseq_ker_nOTU) #6 species for Keratosa

#Haplosclerida
cumporphyseq_hap_nOTU=subset_taxa(cumporphyseq, Order=="Haplosclerida")
#Prune OTUs that are not present in any of the samples
umporphyseq_hap_nOTU=prune_taxa(taxa_sums(cumporphyseq_hap_nOTU) > 0, cumporphyseq_hap_nOTU)
ntaxa(umporphyseq_hap_nOTU) #19 species for Haplosclerida

#Tetractinellida      
cumporphyseq_tet_nOTU=subset_taxa(cumporphyseq, Order=="Tetractinellida")
#Prune OTUs that are not present in any of the samples
umporphyseq_tet_nOTU=prune_taxa(taxa_sums(cumporphyseq_tet_nOTU) > 0, cumporphyseq_tet_nOTU)
ntaxa(umporphyseq_tet_nOTU) #3 species for Tetractinellida

#Suberitida      
cumporphyseq_sub_nOTU=subset_taxa(cumporphyseq, Order=="Suberitida")
#Prune OTUs that are not present in any of the samples
umporphyseq_sub_nOTU=prune_taxa(taxa_sums(cumporphyseq_sub_nOTU) > 0, cumporphyseq_sub_nOTU)
ntaxa(umporphyseq_sub_nOTU) #11 species for Suberitida

#Tethyida     
cumporphyseq_teth_nOTU=subset_taxa(cumporphyseq, Order=="Tethyida")
#Prune OTUs that are not present in any of the samples
umporphyseq_teth_nOTU=prune_taxa(taxa_sums(cumporphyseq_teth_nOTU) > 0, cumporphyseq_teth_nOTU)
ntaxa(umporphyseq_teth_nOTU) #6 species for Tethyida

#Poecilosclerida
cumporphyseq_poec_nOTU=subset_taxa(cumporphyseq, Order=="Poecilosclerida")
#Prune OTUs that are not present in any of the samples
umporphyseq_poec_nOTU=prune_taxa(taxa_sums(cumporphyseq_poec_nOTU) > 0, cumporphyseq_poec_nOTU)
ntaxa(umporphyseq_poec_nOTU) #6 species for Poecilosclerida



#Plotting distribution...
plot_bar(cumporphyseq, x="TMT", fill="Class")
view(cumporphyseq)
#without normalization....
cumbyclass<-psmelt(cumporphyseq) %>%
  ggplot(aes(x = TMT, y = Abundance, fill = Class)) +viridis::scale_color_viridis(discrete = TRUE) +
  theme_minimal() +scale_fill_manual("legend", values=c("Demospongiae" = "#67A9CF","Calcarea" ="#FDBB84","Homoscleromorpha" = "#969696"))+
  theme(legend.position = "bottom")+
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Relative Abundance", title = "Class Relative Abundance")
cumbyclass
head(sample_names(cumporphyseq))

#merge sample by treatment first
cumporphyseqTMTmerge<-merge_samples(cumporphyseq , "TMT")
SD= merge_samples(sample_data(cumporphyseq), "TMT")
print(SD[, "TMT"])
print(cumporphyseqTMTmerge)
sample_names(cumporphyseqTMTmerge)
sample_data(cumporphyseqTMTmerge)$TMT
identical(SD,sample_data(cumporphyseqTMTmerge))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(cumporphyseqTMTmerge)$UNIT <- factor(sample_names(cumporphyseqTMTmerge))
#And change the columns of interest from the metadata to factors as well by doing the following:
sample_data(cumporphyseqTMTmerge)$TMT <- factor(sample_names(cumporphyseqTMTmerge))
#And change the columns of interest from the metadata to factors as well by doing the following:
sample_data(cumporphyseqTMTmerge)$UNIT <- factor(sample_names(cumporphyseqTMTmerge))

####
#heatmaps

```{r heatmap, echo=FALSE}


resorted_taxa<-c("OTU86","OTU84","OTU82","OTU85","OTU83","OTU104","OTU88","OTU87","OTU17","OTU106","OTU21","OTU2","OTU19","OTU16","OTU3","OTU91","OTU14","OTU1","OTU98","OTU15","OTU8","OTU9","OTU23","OTU24","OTU7","OTU10","OTU11","OTU64","OTU65","OTU95","OTU109","OTU67","OTU102","OTU59","OTU60","OTU66","OTU101","OTU62","OTU68","OTU72","OTU71","OTU70","OTU69","OTU74","OTU77","OTU80","OTU73","OTU78","OTU79","OTU204","OTU29","OTU41","OTU33","OTU32","OTU96","OTU31","OTU103","OTU38","OTU40","OTU97","OTU28","OTU27","OTU30","OTU37","OTU26","OTU25","OTU35","OTU111","OTU50","OTU49","OTU57","OTU90","OTU94","OTU56","OTU58","OTU89","OTU54","OTU53","OTU52","OTU55","OTU51","OTU93","OTU46","OTU42","OTU44")

#need to get rid of empty samples
cumporphyseqTMTmerge_pruned=prune_taxa(taxa_sums(cumporphyseqTMTmerge) > 0, cumporphyseqTMTmerge)

cumporphyseqTMTmerge_pruned_hm <- plot_heatmap(cumporphyseqTMTmerge_pruned, "PCoA", "jaccard","TMT", sample.order='TMT', taxa.order= resorted_taxa, taxa.label = "Species", low="#66CCFF", high="#000033", na.value="white") + theme_minimal()
cumporphyseqTMTmerge_pruned_hm
####


#transform for all
cumporphyseqTMTmerge_pruned_normalized_relabun = transform_sample_counts(cumporphyseqTMTmerge_pruned, function(OTU) OTU/sum(OTU) * 100)
#export table with abundance number for species...
cumporphyseqTMTmerge_pruned_normalized_relabun 
cumporphyseqTMTmerge_pruned_normalized_relabun_psmelt<-psmelt(cumporphyseqTMTmerge_pruned_normalized_relabun)
  
head(cumporphyseqTMTmerge_pruned_normalized_relabun_psmelt)

cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Species <-cumporphyseqTMTmerge_pruned_normalized_relabun %>% tax_glom(taxrank = "Species") %>%
  psmelt() %>% select(OTU, Species, TMT, Abundance) %>% 
  spread(TMT, Abundance) 
head(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Species)
#Save data.frame to file
write.table(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Species, file = "cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Species.tsv", sep = "\t", quote = F, 
            row.names = F, col.names = T)
#exporting relaive abundance (%) by class
cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class<-cumporphyseqTMTmerge_pruned_normalized_relabun %>% tax_glom(taxrank = "Class") %>%
  psmelt() %>% select(OTU, Class, TMT, Abundance) %>% 
  spread(TMT, Abundance) 

#Save data.frame to file
write.table(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class, file = "cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class.tsv", sep = "\t", quote = F, 
            row.names = F, col.names = T)


#making a table of %relative abundance of classes
#starting with ATACO2
Class_cum_ATACO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$ATACO2, by=list(Class=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$Class), FUN=sum)
#rename column
colnames(Class_cum_ATACO2_aggregate)[colnames(Class_cum_ATACO2_aggregate) == "x"] <- "Abundance"
view(Class_cum_ATACO2_aggregate)
#adding column (ATACO2)
TMT_ATACO2<-data.frame(TMT=c('ATACO2','ATACO2','ATACO2'))
View(TMT_ATACO2)
view(Class_cum_ATACO2_aggregate)
Class_cum_ATACO2_aggregate<-cbind(TMT_ATACO2,Class_cum_ATACO2_aggregate)
Class_cum_ATACO2_aggregate

#ATHCO2
Class_cum_ATHCO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$ATHCO2, by=list(Class=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$Class), FUN=sum)
#rename column
colnames(Class_cum_ATHCO2_aggregate)[colnames(Class_cum_ATHCO2_aggregate) == "x"] <- "Abundance"
view(Class_cum_ATHCO2_aggregate)
#adding column (ATHCO2)
TMT_ATHCO2<-data.frame(TMT=c('ATHCO2','ATHCO2','ATHCO2'))
View(TMT_ATHCO2)
view(Class_cum_ATHCO2_aggregate)
Class_cum_ATHCO2_aggregate<-cbind(TMT_ATHCO2,Class_cum_ATHCO2_aggregate)
Class_cum_ATHCO2_aggregate

#HTACO2
Class_cum_HTACO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$HTACO2, by=list(Class=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$Class), FUN=sum)
#rename column
colnames(Class_cum_HTACO2_aggregate)[colnames(Class_cum_HTACO2_aggregate) == "x"] <- "Abundance"
view(Class_cum_HTACO2_aggregate)
#adding column (HTACO2)
TMT_HTACO2<-data.frame(TMT=c('HTACO2','HTACO2','HTACO2'))
View(TMT_HTACO2)
view(Class_cum_HTACO2_aggregate)
Class_cum_HTACO2_aggregate<-cbind(TMT_HTACO2,Class_cum_HTACO2_aggregate)
view(Class_cum_HTACO2_aggregate)

#HTHCO2
Class_cum_HTHCO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$HTHCO2, by=list(Class=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Class$Class), FUN=sum)
#rename column
colnames(Class_cum_HTHCO2_aggregate)[colnames(Class_cum_HTHCO2_aggregate) == "x"] <- "Abundance"
view(Class_cum_HTHCO2_aggregate)
#adding column (HTHCO2)
TMT_HTHCO2<-data.frame(TMT=c('HTHCO2','HTHCO2','HTHCO2'))
View(TMT_HTHCO2)
view(Class_cum_HTHCO2_aggregate)
Class_cum_HTHCO2_aggregate<-cbind(TMT_HTHCO2,Class_cum_HTHCO2_aggregate)
view (Class_cum_HTHCO2_aggregate)


Class_cum_allTMTs_aggregate <- Class_cum_HTHCO2_aggregate %>%
  bind_rows(Class_cum_ATHCO2_aggregate) %>%
  bind_rows(Class_cum_ATACO2_aggregate)%>% bind_rows(Class_cum_HTACO2_aggregate)
view(Class_cum_allTMTs_aggregate)
write.csv(Class_cum_allTMTs_aggregate, "Class_cum_allTMTs_aggregate_percentrelabun.csv", row.names = F)

#Class_cum_allTMTs_aggregate_percentrelabun.csv ahs the total abudance percentage for each class across treatments!

#make a table 
plot<-psmelt(cumporphyseqTMTmerge_pruned_normalized_relabun) %>%
  ggplot(aes(x = TMT, y = Abundance, fill = Class)) +viridis::scale_color_viridis(discrete = TRUE)+
  theme_minimal() +scale_fill_manual("legend", values=c("Demospongiae" = "#67A9CF","Calcarea" ="#FDBB84","Homoscleromorpha" = "#969696"))+
  theme(legend.position = "bottom")+
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+  labs(y = "Relative Abundance", title = "Class Relative Abundance")
plot

#exporting relaive abundance (%) by Order


#exporting relaive abundance (%) by Order
cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order<-cumporphyseqTMTmerge_pruned_normalized_relabun %>% tax_glom(taxrank = "Order") %>%
  psmelt() %>% select(OTU, Order, TMT, Abundance) %>% 
  spread(TMT, Abundance) 

view(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order)

#Save data.frame to file
write.table(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order, file = "cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order.tsv", sep = "\t", quote = F, 
            row.names = F, col.names = T)


view(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order)
#making a table of %relative abundance of Order     
#starting with ATACO2
Order_cum_ATACO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$ATACO2, by=list(Order=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$Order), FUN=sum)
#rename column
colnames(Order_cum_ATACO2_aggregate)[colnames(Order_cum_ATACO2_aggregate) == "x"] <- "Abundance"
view(Order_cum_ATACO2_aggregate)
#adding column (ATACO2)
TMT_ATACO2<-data.frame(TMT=c('ATACO2','ATACO2','ATACO2', 'ATACO2','ATACO2','ATACO2','ATACO2','ATACO2','ATACO2','ATACO2','ATACO2'))
View(TMT_ATACO2)
view(Order_cum_ATACO2_aggregate)
Order_cum_ATACO2_aggregate<-cbind(TMT_ATACO2,Order_cum_ATACO2_aggregate)
Order_cum_ATACO2_aggregate

#ATHCO2
Order_cum_ATHCO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$ATHCO2, by=list(Order=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$Order), FUN=sum)
#rename column
colnames(Order_cum_ATHCO2_aggregate)[colnames(Order_cum_ATHCO2_aggregate) == "x"] <- "Abundance"
view(Order_cum_ATHCO2_aggregate)
#adding column (ATHCO2)
TMT_ATHCO2<-data.frame(TMT=c('ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2','ATHCO2'))
View(TMT_ATHCO2)
view(Order_cum_ATHCO2_aggregate)
Order_cum_ATHCO2_aggregate<-cbind(TMT_ATHCO2,Order_cum_ATHCO2_aggregate)
Order_cum_ATHCO2_aggregate

#HTACO2
Order_cum_HTACO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$HTACO2, by=list(Order=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$Order), FUN=sum)
#rename column
colnames(Order_cum_HTACO2_aggregate)[colnames(Order_cum_HTACO2_aggregate) == "x"] <- "Abundance"
view(Order_cum_HTACO2_aggregate)
#adding column (HTACO2)
TMT_HTACO2<-data.frame(TMT=c('HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2','HTACO2'))
View(TMT_HTACO2)
view(Order_cum_HTACO2_aggregate)
Order_cum_HTACO2_aggregate<-cbind(TMT_HTACO2,Order_cum_HTACO2_aggregate)
view(Order_cum_HTACO2_aggregate)

#HTHCO2
Order_cum_HTHCO2_aggregate<-aggregate(cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$HTHCO2, by=list(Order=cumporphyseqTMTmerge_pruned_normalized_relabun_fortablexport_Order$Order), FUN=sum)
#rename column
colnames(Order_cum_HTHCO2_aggregate)[colnames(Order_cum_HTHCO2_aggregate) == "x"] <- "Abundance"
view(Order_cum_HTHCO2_aggregate)
#adding column (HTHCO2)
TMT_HTHCO2<-data.frame(TMT=c('HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2','HTHCO2'))
View(TMT_HTHCO2)
view(Order_cum_HTHCO2_aggregate)
Order_cum_HTHCO2_aggregate<-cbind(TMT_HTHCO2,Order_cum_HTHCO2_aggregate)
view (Order_cum_HTHCO2_aggregate)


Order_cum_allTMTs_aggregate <- Order_cum_HTHCO2_aggregate %>%
  bind_rows(Order_cum_ATHCO2_aggregate) %>%
  bind_rows(Order_cum_ATACO2_aggregate)%>% bind_rows(Order_cum_HTACO2_aggregate)
view(Order_cum_allTMTs_aggregate)
write.csv(Order_cum_allTMTs_aggregate, "Order_cum_allTMTs_aggregate_percentrelabun.csv", row.names = F)


#quantifying abudnacne proportions relative to the control
#subset cumporphyseq for Calcarea only...
cumporphyseq_Calcarea<-subset_taxa(cumporphyseq, Class=="Calcarea")

#subset for Homoscleromorpha only...
cumporphyseq_Homoscleromorpha<-subset_taxa(cumporphyseq, Class=="Homoscleromorpha")

#Subset Keratosa
cumporphyseq_Ker<-subset_taxa(cumporphyseq, Subclass=="Keratosa")

#transform for Haplosclerida
cumporphyseq_hap<-subset_taxa(cumporphyseq, Order=="Haplosclerida")

#transform for Poecilosclerida  
cumporphyseq_poe<-subset_taxa(cumporphyseq, Order=="Poecilosclerida")

#transform for Suberitida  
cumporphyseq_sub<-subset_taxa(cumporphyseq, Order=="Suberitida")

#Tethyida
cumporphyseq_teth<-subset_taxa(cumporphyseq, Order=="Tethyida")

#Tetractinellida
cumporphyseq_tetr<-subset_taxa(cumporphyseq, Order=="Tetractinellida")

#subset for Demospongiae only...
cumporphyseq_Demospongiae<-subset_taxa(cumporphyseq, Class=="Demospongiae")

#transform for Demospongiae
normalized_cumporphyseq_Demospongiae_relabun = transform_sample_counts(cumporphyseq_Demospongiae, function(OTU) OTU/sum(OTU) * 100/6)


#if using Order we should pick the most common Orders.,. so lets do that...
normalized_cumporphyseq_Demospongiae_relabun
normalized_cumporphyseq_Demospongiae_relabun_with_99per<-phyloseq_filter_taxa_tot_fraction(normalized_cumporphyseq_Demospongiae_relabun, frac = 0.0001) #use this function to get rid of a certain fraction

# with order
pordercum<-psmelt(normalized_cumporphyseq_Demospongiae_relabun_with_99per) %>%
  ggplot(aes(x = TMT, y = Abundance, fill = Order))  +theme_minimal() +
  theme(legend.position = "bottom")+
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +scale_fill_manual("legend", values=c("Haplosclerida" = "#FFB053","Tetractinellida" ="#B69C8A","Suberitida" = "#BA9DC6", "Tethyida" = "#E59BA0",
                                                                                                              "Chondrillida" = "#506598","Poecilosclerida" = "#91BC87", "Dendroceratida" = "#8AB8D0", "Dictyoceratida" = "#3F4812"))
  labs(y = "Relative Abundance", title = "Order Relative Abundance within Demospongiae")
  pordercum
#creating sample_sums table for all subgroups to make a proportion comparison figure with the control. 
 #Demospongiae
  sample_sums(cumporphyseq_Demospongiae)
  cumporphyseq_Demospongiae_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_Demospongiae)) 
  write.csv(cumporphyseq_Demospongiae_sumabun_df, "cumporphyseq_Demospongiae_sumabun_df.csv", row.names = F)

#Calcarea
  sample_sums(cumporphyseq_Calcarea)
  cumporphyseq_Calcarea_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_Calcarea)) 
  write.csv(cumporphyseq_Calcarea_sumabun_df, "cumporphyseq_Calcarea_sumabun_df.csv", row.names = F)

#Homoscleromorpha
  sample_sums(cumporphyseq_Homoscleromorpha)
  cumporphyseq_Homo_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_Homoscleromorpha)) 
  write.csv(cumporphyseq_Homo_sumabun_df, "cumporphyseq_Homo_sumabun_df.csv", row.names = F)
  
#Keratosa
  sample_sums(cumporphyseq_Ker)
  cumporphyseq_Ker_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_Ker)) 
  write.csv(cumporphyseq_Ker_sumabun_df, "cumporphyseq_Ker_sumabun_df.csv", row.names = F)

#Haploclerida 
  sample_sums(cumporphyseq_hap)
  cumporphyseq_hap_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_hap)) 
  View(cumporphyseq_hap_sumabun_df)
  write.csv(cumporphyseq_hap_sumabun_df, "cumporphyseq_hap_sumabun_df.csv", row.names = F)
  View(cumporphyseq_hap)
#Poecilosclerida
  sample_sums(cumporphyseq_poe)
  cumporphyseq_poe_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_poe)) 
  write.csv(cumporphyseq_poe_sumabun_df, "cumporphyseq_poe_sumabun_df.csv", row.names = F)
  
#Suberitida
  sample_sums(cumporphyseq_sub)
  cumporphyseq_sub_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_sub)) 
  View(cumporphyseq_sub_sumabun_df)
  write.csv(cumporphyseq_sub_sumabun_df, "cumporphyseq_sub_sumabun_df.csv", row.names = F)
  
#Tethyida
  sample_sums(cumporphyseq_teth)
  cumporphyseq_teth_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_teth)) 
  write.csv(cumporphyseq_teth_sumabun_df, "cumporphyseq_teth_sumabun_df.csv", row.names = F)
  
#Tetractinellida
  sample_sums(cumporphyseq_tetr)
  cumporphyseq_tetr_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_tetr)) 
  write.csv(cumporphyseq_tetr_sumabun_df, "cumporphyseq_tetr_sumabun_df.csv", row.names = F)
  
#calculated proportion differences in excel and now importing into R to plot. 
  cumprop<-read.csv("cumprop.csv") #import .csv
  cumprop[-c(28), ] #delete fow 28

#plot dots for proportions
  
  ratiop2 <- cumprop %>%
    mutate%>%ggplot(aes(Ratio, Sponge.group, colour = factor(TMT)))+geom_point(size = 4)+theme_classic()+ scale_color_manual(values = c("HTACO2" ="#D55E00","ATHCO2" = "#56B4E9", "HTHCO2" = "#999999"))+scale_y_discrete(name ="Sopnge.group",limits=c("Haplosclerida", "Poecilosclerida", "Suberitida", "Tethyida", "Tetractinellida", "Keratosa", "Homoscleromorpha", "Calcarea", "Demospongiae"))
ratiop2

######
#venn

