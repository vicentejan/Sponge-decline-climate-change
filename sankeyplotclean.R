#making a sankey diagram....

#calculating diversity of Porifera based on mergingUNITTime. The diversity of each unit comprises that of every plate and side. 
```{r setup, include=FALSE}
#merge samples by unitTime. This way we can analyze diversity of entire Unit at each time point. 

physeqpor=physeq

#Prune OTUs that are not present in any of the samples
physeqpor=prune_taxa(taxa_sums(physeq) > 0, physeq) #This actually gets rid of samples that have more than 0 so not exactly sure what prune_taxa is doing. Need to figure this out... Made a heatmap with and without


mergedTMTTime = merge_samples(physeq, "TMTTime")
SD= merge_samples(sample_data(physeq), "TMTTime")
print(SD[, "TMTTime"])
print(mergedTMTTime)
sample_names(mergedTMTTime)
sample_data(mergedTMTTime)$TMTTime
identical(SD,sample_data(mergedTMTTime))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(mergedTMTTime)$UNITTime <- factor(sample_names(mergedTMTTime))
#And change the columns of interest from the metadata to factors as well by doing the following:
sample_data(mergedTMTTime)$UNITTime <- factor(sample_names(mergedTMTTime))
#And change the columns of interest from the metadata to factors as well by doing the following:
sample_data(mergedTMTTime)$TMTTime <- factor(sample_names(mergedTMTTime))
#And change the columns of interest from the metadata to factors as well by doing the following:
#Changing the TMT values
sample_data(mergedTMTTime)$TMT<-factor(sample_data(mergedTMTTime)$TMT)
#Then Manually inserting the original names of the factors by doing the following:
sample_data(mergedTMTTime)$Date<-factor(sample_data(mergedTMTTime)$Date)
sample_data(mergedTMTTime)$TMT <- revalue(sample_data(mergedTMTTime)$TMT, c("1"="ATACO2", "2"="ATHCO2", "3"="HTACO2", "4"="HTHCO2", "5"="Intake"))

#Changing the Date values to factor
sample_data(mergedTMTTime)$Date <- factor(sample_names(mergedTMTTime))
#Didn't change time.
mergedDATE<-read_excel("DATE_TMTTime.xlsx", 1) #made a column in excel with the correct date
sample_data(mergedTMTTime)$DATE<-mergedDATE$TIME #added date column 
sample_data(mergedTMTTime)$DATE <- factor(sample_data(mergedTMTTime)$DATE) #Changing the Date values to factor
mergedTMTTime_nointake=subset_samples(mergedTMTTime, TMT != "Intake")#now that I have mergedUNITTime I can eliminate intake and work on treatments only
View(mergedTMTTime_nointake)

mergedTMTTime_nointake_por<-subset_taxa(mergedTMTTime_nointake, mergedTMTTime_nointake@tax_table[, 1]=="Porifera")


#sankey plot #downloaded install_github("adrientaudiere/Miscmetabar")
#> Also defined by ‘RNeXML’
skplotOTU<-if (requireNamespace("networkD3")) {
  sankey_pq(mergedTMTTime_nointake_por, fact = "TMT")}
#this plots the number of OTUs but 
skplotOTU  #This lists the OTUs

saveNetwork(skplotOTU, "skplotOTU.html") #saves plot in html format

webshot::webshot("skplotOTU.html","skplotOTU.png",          vwidth = 1000, vheight = 900, zoom = 2)
webshot::webshot("skplotOTU.html","skplotOTU.pdf",          vwidth = 1000, vheight = 900, zoom = 2)

#but with add_nb_seq = TRUE it adds the number of observations for each OTU and.  
skplotOTUabun<-if (requireNamespace("networkD3")) {
  sankey_pq(cumporphyseq, fact = "TMT", add_nb_seq = TRUE)
}
skplotOTUabun
saveNetwork(skplotOTUabun, "skplotOTUabun.html") #saves plot in html format
webshot::webshot("skplotOTUabun.html","skplotOTUabun.pdf",          vwidth = 1000, vheight = 900, zoom = 2)



#issue solved by using Webshhot2
ggsave()
webshot::webshot("skplotOTUabun.html","skplotOTUabun.svg",          vwidth = 1000, vheight = 900, zoom = 2)
webshot("skplotOTUabun.svg", vwidth = 938, vheight = 821)

webshot::webshot("skplotOTUabun.html","skplotOTUabun.svg", vwidth = 938, vheight = 821)


#double checking abundance values....
#transform for Suberitida  
cumporphyseq_sub_sandkey<-subset_taxa(mergedTMTTime_nointake_por, Order=="Suberitida")
#Suberitida
sample_sums(cumporphyseq_sub_sandkey)
cumporphyseq_sub_sandkey_sub_sumabun_df<-as.data.frame(sample_sums(cumporphyseq_sub_sandkey)) 
View(cumporphyseq_sub_sandkey_sub_sumabun_df)
write.csv(cumporphyseq_sub_sumabun_df, "cumporphyseq_sub_sumabun_df.csv", row.names = F)
### there is clearly a mistake so lets try it again using the dataset already gathered from
#the heatmamp and cumulative porifera...

#Haplosclerida
cumporphyseq_hap_sandkey<-subset_taxa(mergedTMTTime_nointake_por, Order=="Haplosclerida")
sample_sums(cumporphyseq_hap_sandkey)
cumporphyseq_hap_sumabun_sk_df<-as.data.frame(sample_sums(cumporphyseq_hap_sandkey)) 
View(cumporphyseq_hap_sumabun_sk_df)
write.csv(cumporphyseq_sub_sumabun_df, "cumporphyseq_sub_sumabun_df.csv", row.names = F)

#Double  checked data by running (cumporphyseqs) the Sankey plot again with the phyloseq object from the cumulative diversity section and the same plot was produced
#there is no error in the data. 

Doing a Sankey with "cumporphyseqs" which was the phyloseq object created for cumulative abundance. 

skplotC<-if (requireNamespace("networkD3")) {
  sankey_pq(cumporphyseq, fact = "TMT")}

skplotC
saveNetwork(skplot, "skp.html") #saves plot in html format

