#Make a richness output for diversity so that we can plot the different metrics. 
#Porifera

#calculating diversity of Porifera based on mergingUNITTime. The diversity of each unit comprises that of every plate and side. 
```{r setup, include=FALSE}
#merge samples by unitTime. This way we can analyze diversity of entire Unit at each time point. 

physeqpor=physeq

#Prune OTUs that are not present in any of the samples
physeqpor=prune_taxa(taxa_sums(physeq) > 0, physeq) #This actually gets rid of samples that have more than 0 so not exactly sure what prune_taxa is doing. Need to figure this out... Made a heatmap with and without

mergedUNITTime = merge_samples(physeqpor, "UNITTime")
SD= merge_samples(sample_data(physeqpor), "UNITTime")
print(SD[, "UNITTime"]) #Error in validObject(.Object) : 
#invalid class “sample_data” object: Sample Data must have non-zero dimensions.
print(mergedUNITTime)
sample_names(mergedUNITTime)
sample_data(mergedUNITTime)$Date
identical(SD,sample_data(mergedUNITTime))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(mergedUNITTime)$UNITTime <- factor(sample_names(mergedUNITTime))
#And change the columns of interest from the metadata to factors as well by doing the following:
#Changing the TMT values
sample_data(mergedUNITTime)$TMT<-factor(sample_data(mergedUNITTime)$TMT)
#Then Manually inserting the original names of the factors by doing the following:
sample_data(mergedUNITTime)$Date<-factor(sample_data(mergedUNITTime)$Date)
sample_data(mergedUNITTime)$TMT <- revalue(sample_data(mergedUNITTime)$TMT, c("1"="ATACO2", "2"="ATHCO2", "3"="HTACO2", "4"="HTHCO2", "5"="Intake"))
#Changing the Date values to factor
sample_data(mergedUNITTime)$Date <- factor(sample_names(mergedUNITTime))
#Didn't change time.

mergedDATE<-read_excel("DATE.xlsx", 1) #made a column in excell with the correct date
sample_data(mergedUNITTime)$DATE<-mergedDATE$TIME #added date column 
sample_data(mergedUNITTime)$DATE <- factor(sample_data(mergedUNITTime)$DATE) #Changing the Date values to factor
sample_data(mergedUNITTime)$pco2<-sample_data(mergedUNITTime)$TMT #added pco2 column 
sample_data(mergedUNITTime)$pco2 <- revalue(sample_data(mergedUNITTime)$pco2, c("ATACO2"="L", "ATHCO2"="H", "HTACO2"="L", "HTHCO2"="H", "Intake"="L"))

sample_data(mergedUNITTime)$temp<-sample_data(mergedUNITTime)$TMT #added temperature column 
sample_data(mergedUNITTime)$temp <- revalue(sample_data(mergedUNITTime)$temp, c("ATACO2"="L", "ATHCO2"="L", "HTACO2"="H", "HTHCO2"="H", "Intake"="L"))

#subset to last time point
Last_mergedUNITTime=subset_samples(mergedUNITTime, DATE == "2018-06-11")
Calc_mergedUNITTime=subset_taxa(mergedUNITTime, Class == "Calcarea")
```

#Rarefraction by TMT=ATACO2



#Rarefraction by date "Date". Samples from "ATACO2", HTACO2", "HTHCO2", "ATHC02" separately subseted from phyloseq


```{r setup, include=FALSE}
#Subset data by TMT="ATACO2" 
#Subset data by TMT="Ambient" 
Ambientphyseqpor=subset_samples(physeqpor,TMT=="ATACO2")
AmbientphyseqpormergedDate = merge_samples(Ambientphyseqpor, "Date")
SD= merge_samples(sample_data(Ambientphyseqpor), "Date")
print(SD[, "Date"])
print(AmbientphyseqpormergedDate)
sample_names(AmbientphyseqpormergedDate)
sample_data(AmbientphyseqpormergedDate)$Date
identical(SD,sample_data(AmbientphyseqpormergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqpormergedDate)$Date <- factor(sample_names(AmbientphyseqpormergedDate))
View(AmbientphyseqpormergedDate)
#Getting abundance
poriferaAmbient_abund_mergedDate=as(otu_table(AmbientphyseqpormergedDate), "matrix")
View(poriferaAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poriferaAmbient_abund_mergedDate_abs_pre<-poriferaAmbient_abund_mergedDate
write.csv(poriferaAmbient_abund_mergedDate, "poriferaAmbient_abund_mergedDate.csv", row.names = F)
poriferaAmbient_abund_mergedDate_abs_pre[poriferaAmbient_abund_mergedDate_abs_pre >0] <-1 #making a presence absence data frame
View(poriferaAmbient_abund_mergedDate_abs_pre)
write.csv(poriferaAmbient_abund_mergedDate_abs_pre, "poriferaAmbient_abund_mergedDate_abs_pre.csv", row.names = F)

###subsetting samples by "ambient" and then merging sample by "Unit"" 
#which will combine the total diversity of the experiment for this treatment and calculate mean values of diversity#####
Ambientphyseqpormergedunit = merge_samples(Ambientphyseqpor, "UNIT")
SD= merge_samples(sample_data(Ambientphyseqpor), "UNIT")
print(SD[, "UNIT"])
print(Ambientphyseqpormergedunit)
sample_names(Ambientphyseqpormergedunit)
sample_data(Ambientphyseqpormergedunit)$UNIT
identical(SD,sample_data(Ambientphyseqpormergedunit))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(Ambientphyseqpormergedunit)$UNIT <- factor(sample_names(Ambientphyseqpormergedunit))
#Determine diversity within Porifera Ambient.
Richness_Ambientphyseqpormergedunit<-estimate_richness(Ambientphyseqpormergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqpormergedunit, "Richness_Ambientphyseqpormergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqpormergedunit)=sample_names(Ambientphyseqpormergedunit)
Richness_Ambientphyseqpormergedunit<-dfRowName(Richness_Ambientphyseqpormergedunit, name ="UNIT")
View(Richness_Ambientphyseqpormergedunit)
#need to add a column "TMT" to indicate "Ambient"
Richness_Ambientphyseqpormergedunit<-cbind(Richness_Ambientphyseqpormergedunit, TMT='ATACO2')
View(Richness_Ambientphyseqpormergedunit)
#Calculating mean values of Observed, SHannon, Simpson for the Ambient. which will be shown in table xx diversity. 
#Observed
Richness_Ambientphyseqpormergedunit_sum<-summarySE(Richness_Ambientphyseqpormergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqpormergedunit_sum)
#Shannon
Shannon_Ambientphyseqpormergedunit_sum<-summarySE(Richness_Ambientphyseqpormergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqpormergedunit_sum)
#Simpson
Simpson_Ambientphyseqpormergedunit_sum<-summarySE(Richness_Ambientphyseqpormergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqpormergedunit_sum)
#Getting abundance
poriferaAmbient_abund_mergedunit=as(otu_table(Ambientphyseqpormergedunit), "matrix")
View(poriferaAmbient_abund_mergedunit)
head(poriferaAmbient_abund_mergedunit)
write.csv(poriferaAmbient_abund_mergedunit, "poriferaAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
poriferaAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(poriferaAmbient_abund_mergedunit), data.frame(poriferaAmbient_abund_mergedunit, row.names=NULL))
View(poriferaAmbient_abund_mergedunit_rn_to_col)
head(poriferaAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poriferaAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
poriferaAmbient_abund_mergedunit_rn_to_col_asdf<-as_data_frame(poriferaAmbient_abund_mergedunit) 
View(poriferaAmbient_abund_mergedunit_rn_to_col_asdf)
names(poriferaAmbient_abund_mergedunit_rn_to_col_asdf)
#Sum rows
poriferaAmbient_abund_mergedunit_rn_to_col_total<-cbind(poriferaAmbient_abund_mergedunit_rn_to_col_asdf, total =rowSums(poriferaAmbient_abund_mergedunit_rn_to_col_asdf)) 
View(poriferaAmbient_abund_mergedunit_rn_to_col_total)
names(poriferaAmbient_abund_mergedunit_rn_to_col_total)
#add column "Unit" from "poriferaAmbient_abund_mergedunit_rn_to_col" to "poriferaAmbient_abund_mergedunit_rn_to_col_total"
poriferaAmbient_abund_mergedunit_rn_to_col_total$Unit<-poriferaAmbient_abund_mergedunit_rn_to_col$Unit
View(poriferaAmbient_abund_mergedunit_rn_to_col_total)
names(poriferaAmbient_abund_mergedunit_rn_to_col_total)
#Eliminate column 1 through 94. 
poriferaAmbient_abund_mergedunit_rn_to_col_total_clean<-poriferaAmbient_abund_mergedunit_rn_to_col_total[,-(1:93)]
View(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean)
head(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean)
#need to add a column "TMT" to indicate "Ambient"
poriferaAmbient_abund_mergedunit_rn_to_col_total_clean<-cbind(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean, TMT='ATACO2')
View(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean)
head(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean)
names(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean)
write.csv(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean, "poriferaAmbient_abund_mergedunit_rn_to_col_total_clean.csv", row.names = F)
#SUmmarySE
poriferaAmbient_abund_mergedunit_rn_to_col_total_clean_sum<-summarySE(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poriferaAmbient_abund_mergedunit_rn_to_col_total_clean_sum)

#Now calculate diversity within the different sponge subgroups within "Ambientphyseqpormergedunit"" starting with Class Homoscleromorpha
Ambientphyseqhomomergedunit<-subset_taxa(Ambientphyseqpormergedunit, Class=="Homoscleromorpha")
Richness_Ambientphyseqhomomergedunit<-estimate_richness(Ambientphyseqhomomergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqhomomergedunit, "Richness_Ambientphyseqhomomergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqhomomergedunit)=sample_names(Ambientphyseqhomomergedunit)
Richness_Ambientphyseqhomomergedunit<-dfRowName(Richness_Ambientphyseqhomomergedunit, name ="UNIT")
View(Richness_Ambientphyseqhomomergedunit)
#Observed
Richness_Ambientphyseqhomomergedunit<-cbind(Richness_Ambientphyseqhomomergedunit, TMT='ATACO2')
Richness_Ambientphyseqhomomergedunit_sum<-summarySE(Richness_Ambientphyseqhomomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqhomomergedunit_sum)
#Shannon
Shannon_Ambientphyseqhomomergedunit_sum<-summarySE(Richness_Ambientphyseqhomomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqhomomergedunit_sum)
#Simpson
Simpson_Ambientphyseqhomomergedunit_sum<-summarySE(Richness_Ambientphyseqhomomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqhomomergedunit_sum)
#Getting abundance
HomoAmbient_abund_mergedunit=as(otu_table(Ambientphyseqhomomergedunit), "matrix")
View(HomoAmbient_abund_mergedunit)
write.csv(HomoAmbient_abund_mergedunit, "HomoAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
HomoAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(HomoAmbient_abund_mergedunit), data.frame(HomoAmbient_abund_mergedunit, row.names=NULL))
View(HomoAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(HomoAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
HomoAmbient_abund_mergedunit_asdf<-as_data_frame(HomoAmbient_abund_mergedunit) 
View(HomoAmbient_abund_mergedunit_asdf)
#Sum rows
HomoAmbient_abund_mergedunit_asdf_total<-cbind(HomoAmbient_abund_mergedunit_asdf, total =rowSums(HomoAmbient_abund_mergedunit_asdf)) 
View(HomoAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "HomoAmbient_abund_mergedunit_rn_to_col" to "HomoAmbient_abund_mergedunit_rn_to_col_total"
HomoAmbient_abund_mergedunit_asdf_total$Unit<-HomoAmbient_abund_mergedunit_rn_to_col$Unit
View(HomoAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 8. 
HomoAmbient_abund_mergedunit_rn_to_col_total_clean<-HomoAmbient_abund_mergedunit_asdf_total[,-(1:8)]
View(HomoAmbient_abund_mergedunit_rn_to_col_total_clean)
#need to add a column "TMT" to indicate "Ambient"
HomoAmbient_abund_mergedunit_rn_to_col_total_clean<-cbind(HomoAmbient_abund_mergedunit_rn_to_col_total_clean, TMT='ATACO2')
View(HomoAmbient_abund_mergedunit_rn_to_col_total_clean)
write.csv(HomoAmbient_abund_mergedunit_rn_to_col_total_clean, "HomoAmbient_abund_mergedunit_rn_to_col_total_clean.csv", row.names = F)
#SUmmarySE
HomoAmbient_abund_mergedunit_rn_to_col_total_clean_sum<-summarySE(HomoAmbient_abund_mergedunit_rn_to_col_total_clean, measurevar = "total", groupvars=c("TMT"))
View(HomoAmbient_abund_mergedunit_rn_to_col_total_clean_sum)
head(HomoAmbient_abund_mergedunit_rn_to_col_total_clean_sum)

Class==Demospongiae
Ambientphyseqdemmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Class=="Demospongiae")
Richness_Ambientphyseqdemmergedunit<-estimate_richness(Ambientphyseqdemmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqdemmergedunit, "Richness_Ambientphyseqdemmergedunit.csv", 
          row.names = F)
row.names(Richness_Ambientphyseqdemmergedunit)=sample_names(Ambientphyseqdemmergedunit)
Richness_Ambientphyseqdemmergedunit<-dfRowName(Richness_Ambientphyseqdemmergedunit, name ="UNIT")
View(Richness_Ambientphyseqdemmergedunit)
#Observed
Richness_Ambientphyseqdemmergedunit<-cbind(Richness_Ambientphyseqdemmergedunit, TMT='ATACO2')
Richness_Ambientphyseqdemmergedunit_sum<-summarySE(Richness_Ambientphyseqdemmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqdemmergedunit_sum)
#Shannon
Shannon_Ambientphyseqdemmergedunit_sum<-summarySE(Richness_Ambientphyseqdemmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqdemmergedunit_sum)
#Simpson
Simpson_Ambientphyseqdemmergedunit_sum<-summarySE(Richness_Ambientphyseqdemmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqdemmergedunit_sum)
#Getting abundance
demAmbient_abund_mergedunit=as(otu_table(Ambientphyseqdemmergedunit), "matrix")
names(demAmbient_abund_mergedunit)
View(demAmbient_abund_mergedunit)
write.csv(demAmbient_abund_mergedunit, "demAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
demAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(demAmbient_abund_mergedunit), data.frame(demAmbient_abund_mergedunit, row.names=NULL))
View(demAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(demAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
demAmbient_abund_mergedunit_asdf<-as_data_frame(demAmbient_abund_mergedunit) 
View(demAmbient_abund_mergedunit_asdf)
#Sum rows
demAmbient_abund_mergedunit_asdf_total<-cbind(demAmbient_abund_mergedunit_asdf, total =rowSums(demAmbient_abund_mergedunit_asdf)) 
View(demAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "demAmbient_abund_mergedunit_rn_to_col" to "demAmbient_abund_mergedunit_rn_to_col_total"
demAmbient_abund_mergedunit_asdf_total$Unit<-demAmbient_abund_mergedunit_rn_to_col$Unit
View(demAmbient_abund_mergedunit_asdf_total)
names(demAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 67. 
demAmbient_abund_mergedunit_rn_to_col_total_clean<-demAmbient_abund_mergedunit_asdf_total[,-(1:66)]
View(demAmbient_abund_mergedunit_rn_to_col_total_clean)
#need to add a column "TMT" to indicate "Ambient"
demAmbient_abund_mergedunit_rn_to_col_total_clean<-cbind(demAmbient_abund_mergedunit_rn_to_col_total_clean, TMT='ATACO2')
View(demAmbient_abund_mergedunit_rn_to_col_total_clean)
write.csv(demAmbient_abund_mergedunit_rn_to_col_total_clean, "demAmbient_abund_mergedunit_rn_to_col_total_clean.csv", row.names = F)
#SUmmarySE
demAmbient_abund_mergedunit_rn_to_col_total_clean_sum<-summarySE(demAmbient_abund_mergedunit_rn_to_col_total_clean, measurevar = "total", groupvars=c("TMT"))
View(demAmbient_abund_mergedunit_rn_to_col_total_clean_sum)

#Class==Calcarea
Ambientphyseqcalcmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Class=="Calcarea")
Richness_Ambientphyseqcalcmergedunit<-estimate_richness(Ambientphyseqcalcmergedunit, measures = c("Observed", "Shannon", "Simpson"))
row.names(Richness_Ambientphyseqcalcmergedunit)=sample_names(Ambientphyseqcalcmergedunit)
Richness_Ambientphyseqcalcmergedunit<-dfRowName(Richness_Ambientphyseqcalcmergedunit, name ="UNIT")
View(Richness_Ambientphyseqcalcmergedunit)
Richness_Ambientphyseqcalcmergedunit<-cbind(Richness_Ambientphyseqcalcmergedunit, TMT='ATACO2')
#export the replicate values for richness for t test analysis.
write.csv(Richness_Ambientphyseqcalcmergedunit, "Richness_Ambientphyseqcalcmergedunit.csv", row.names = F)
#Observed
Richness_Ambientphyseqcalcmergedunit_sum<-summarySE(Richness_Ambientphyseqcalcmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqcalcmergedunit_sum)
#Shannon
Shannon_Ambientphyseqcalcmergedunit_sum<-summarySE(Richness_Ambientphyseqcalcmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqcalcmergedunit_sum)
#Simpson
Simpson_Ambientphyseqcalcmergedunit_sum<-summarySE(Richness_Ambientphyseqcalcmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqcalcmergedunit_sum)
#Getting abundance
CalAmbient_abund_mergedunit=as(otu_table(Ambientphyseqcalcmergedunit), "matrix")
View(CalAmbient_abund_mergedunit)
write.csv(CalAmbient_abund_mergedunit, "CalAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
CalAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(CalAmbient_abund_mergedunit), data.frame(CalAmbient_abund_mergedunit, row.names=NULL))
View(CalAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(CalAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
CalAmbient_abund_mergedunit_asdf<-as_data_frame(CalAmbient_abund_mergedunit) 
View(CalAmbient_abund_mergedunit_asdf)
#Sum rows
CalAmbient_abund_mergedunit_asdf_total<-cbind(CalAmbient_abund_mergedunit_asdf, total =rowSums(CalAmbient_abund_mergedunit_asdf)) 
View(CalAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "CalAmbient_abund_mergedunit_rn_to_col" to "CalAmbient_abund_mergedunit_rn_to_col_total"
CalAmbient_abund_mergedunit_asdf_total$Unit<-CalAmbient_abund_mergedunit_rn_to_col$Unit
View(CalAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 19. 
CalAmbient_abund_mergedunit_asdf_total_clean<-CalAmbient_abund_mergedunit_asdf_total[,-(1:19)]
View(CalAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
CalAmbient_abund_mergedunit_asdf_total_clean<-cbind(CalAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(CalAmbient_abund_mergedunit_asdf_total_clean)
write.csv(CalAmbient_abund_mergedunit_asdf_total_clean, "CalAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
CalAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(CalAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(CalAmbient_abund_mergedunit_asdf_total_clean_sum)

#Subclass==Keratosa
Ambientphyseqkermergedunit<-subset_taxa(Ambientphyseqpormergedunit, Subclass=="Keratosa")
Richness_Ambientphyseqkermergedunit<-estimate_richness(Ambientphyseqkermergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqkermergedunit, "Richness_Ambientphyseqkermergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqkermergedunit)=sample_names(Ambientphyseqkermergedunit)
Richness_Ambientphyseqkermergedunit<-dfRowName(Richness_Ambientphyseqkermergedunit, name ="UNIT")
View(Richness_Ambientphyseqkermergedunit)
Richness_Ambientphyseqkermergedunit<-cbind(Richness_Ambientphyseqkermergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqkermergedunit_sum<-summarySE(Richness_Ambientphyseqkermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqkermergedunit_sum)
#Shannon
Shannon_Ambientphyseqkermergedunit_sum<-summarySE(Richness_Ambientphyseqkermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqkermergedunit_sum)
#Simpson
Simpson_Ambientphyseqkermergedunit_sum<-summarySE(Richness_Ambientphyseqkermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqkermergedunit_sum)
#Getting abundance
kerAmbient_abund_mergedunit=as(otu_table(Ambientphyseqkermergedunit), "matrix")
View(kerAmbient_abund_mergedunit)
#Convert Rowname into first column 
kerAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(kerAmbient_abund_mergedunit), data.frame(kerAmbient_abund_mergedunit, row.names=NULL))
View(kerAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(kerAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
kerAmbient_abund_mergedunit_asdf<-as_data_frame(kerAmbient_abund_mergedunit) 
View(kerAmbient_abund_mergedunit_asdf)
#Sum rows
kerAmbient_abund_mergedunit_asdf_total<-cbind(kerAmbient_abund_mergedunit_asdf, total =rowSums(kerAmbient_abund_mergedunit_asdf)) 
View(kerAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "kerAmbient_abund_mergedunit_rn_to_col" to "kerAmbient_abund_mergedunit_rn_to_col_total"
kerAmbient_abund_mergedunit_asdf_total$Unit<-kerAmbient_abund_mergedunit_rn_to_col$Unit
View(kerAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
kerAmbient_abund_mergedunit_asdf_total_clean<-kerAmbient_abund_mergedunit_asdf_total[,-(1:6)]
View(kerAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
kerAmbient_abund_mergedunit_asdf_total_clean<-cbind(kerAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(kerAmbient_abund_mergedunit_asdf_total_clean)
write.csv(kerAmbient_abund_mergedunit_asdf_total_clean, "kerAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
kerAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(kerAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(kerAmbient_abund_mergedunit_asdf_total_clean_sum)

#Order==Suberitida
Ambientphyseqsubmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Order=="Suberitida")
Richness_Ambientphyseqsubmergedunit<-estimate_richness(Ambientphyseqsubmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqsubmergedunit, "Richness_Ambientphyseqsubmergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqsubmergedunit)=sample_names(Ambientphyseqsubmergedunit)
Richness_Ambientphyseqsubmergedunit<-dfRowName(Richness_Ambientphyseqsubmergedunit, name ="UNIT")
View(Richness_Ambientphyseqsubmergedunit)
Richness_Ambientphyseqsubmergedunit<-cbind(Richness_Ambientphyseqsubmergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqsubmergedunit_sum<-summarySE(Richness_Ambientphyseqsubmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqsubmergedunit_sum)
#Shannon
Shannon_Ambientphyseqsubmergedunit_sum<-summarySE(Richness_Ambientphyseqsubmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqsubmergedunit_sum)
#Simpson
Simpson_Ambientphyseqsubmergedunit_sum<-summarySE(Richness_Ambientphyseqsubmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqsubmergedunit_sum)

#Getting abundance
subAmbient_abund_mergedunit=as(otu_table(Ambientphyseqsubmergedunit), "matrix")
View(subAmbient_abund_mergedunit)
write.csv(subAmbient_abund_mergedunit, "subAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
subAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(subAmbient_abund_mergedunit), data.frame(subAmbient_abund_mergedunit, row.names=NULL))
View(subAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(subAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
subAmbient_abund_mergedunit_asdf<-as_data_frame(subAmbient_abund_mergedunit) 
View(subAmbient_abund_mergedunit_asdf)
#Sum rows
subAmbient_abund_mergedunit_asdf_total<-cbind(subAmbient_abund_mergedunit_asdf, total =rowSums(subAmbient_abund_mergedunit_asdf)) 
View(subAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "subAmbient_abund_mergedunit_rn_to_col" to "subAmbient_abund_mergedunit_rn_to_col_total"
subAmbient_abund_mergedunit_asdf_total$Unit<-subAmbient_abund_mergedunit_rn_to_col$Unit
View(subAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
subAmbient_abund_mergedunit_asdf_total_clean<-subAmbient_abund_mergedunit_asdf_total[,-(1:15)]
View(subAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
subAmbient_abund_mergedunit_asdf_total_clean<-cbind(subAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(subAmbient_abund_mergedunit_asdf_total_clean)
write.csv(subAmbient_abund_mergedunit_asdf_total_clean, "subAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
subAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(subAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(subAmbient_abund_mergedunit_asdf_total_clean_sum)

#Tetractinellida
Ambientphyseqtetrmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Order=="Tetractinellida")
Richness_Ambientphyseqtetrmergedunit<-estimate_richness(Ambientphyseqtetrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqtetrmergedunit, "Richness_Ambientphyseqtetrmergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqtetrmergedunit)=sample_names(Ambientphyseqtetrmergedunit)
Richness_Ambientphyseqtetrmergedunit<-dfRowName(Richness_Ambientphyseqtetrmergedunit, name ="UNIT")
View(Richness_Ambientphyseqtetrmergedunit)
Richness_Ambientphyseqtetrmergedunit<-cbind(Richness_Ambientphyseqtetrmergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqtetrmergedunit_sum<-summarySE(Richness_Ambientphyseqtetrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqtetrmergedunit_sum)
#Shannon
Shannon_Ambientphyseqtetrmergedunit_sum<-summarySE(Richness_Ambientphyseqtetrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqtetrmergedunit_sum)
#Simpson
Simpson_Ambientphyseqtetrmergedunit_sum<-summarySE(Richness_Ambientphyseqtetrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqtetrmergedunit_sum)
#Getting abundance
tetrAmbient_abund_mergedunit=as(otu_table(Ambientphyseqtetrmergedunit), "matrix")
View(tetrAmbient_abund_mergedunit)
write.csv(tetrAmbient_abund_mergedunit, "tetrAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
tetrAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(tetrAmbient_abund_mergedunit), data.frame(tetrAmbient_abund_mergedunit, row.names=NULL))
View(tetrAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tetrAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
tetrAmbient_abund_mergedunit_asdf<-as_data_frame(tetrAmbient_abund_mergedunit) 
View(tetrAmbient_abund_mergedunit_asdf)
#Sum rows
tetrAmbient_abund_mergedunit_asdf_total<-cbind(tetrAmbient_abund_mergedunit_asdf, total =rowSums(tetrAmbient_abund_mergedunit_asdf)) 
View(tetrAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "tetrAmbient_abund_mergedunit_rn_to_col" to "tetrAmbient_abund_mergedunit_rn_to_col_total"
tetrAmbient_abund_mergedunit_asdf_total$Unit<-tetrAmbient_abund_mergedunit_rn_to_col$Unit
View(tetrAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
tetrAmbient_abund_mergedunit_asdf_total_clean<-tetrAmbient_abund_mergedunit_asdf_total[,-(1:3)]
View(tetrAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
tetrAmbient_abund_mergedunit_asdf_total_clean<-cbind(tetrAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(tetrAmbient_abund_mergedunit_asdf_total_clean)
write.csv(tetrAmbient_abund_mergedunit_asdf_total_clean, "tetrAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tetrAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(tetrAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tetrAmbient_abund_mergedunit_asdf_total_clean_sum)

#Order==Poecilosclerida
Ambientphyseqpoecrmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Order=="Poecilosclerida")
Richness_Ambientphyseqpoecrmergedunit<-estimate_richness(Ambientphyseqpoecrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqpoecrmergedunit, "Richness_Ambientphyseqpoecrmergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqpoecrmergedunit)=sample_names(Ambientphyseqpoecrmergedunit)
Richness_Ambientphyseqpoecrmergedunit<-dfRowName(Richness_Ambientphyseqpoecrmergedunit, name ="UNIT")
View(Richness_Ambientphyseqpoecrmergedunit)
Richness_Ambientphyseqpoecrmergedunit<-cbind(Richness_Ambientphyseqpoecrmergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqpoecrmergedunit_sum<-summarySE(Richness_Ambientphyseqpoecrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqpoecrmergedunit_sum)
#Shannon
Shannon_Ambientphyseqpoecrmergedunit_sum<-summarySE(Richness_Ambientphyseqpoecrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqpoecrmergedunit_sum)
#Simpson
Simpson_Ambientphyseqpoecrmergedunit_sum<-summarySE(Richness_Ambientphyseqpoecrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqpoecrmergedunit_sum)
#Getting abundance
poecrAmbient_abund_mergedunit=as(otu_table(Ambientphyseqpoecrmergedunit), "matrix")
View(poecrAmbient_abund_mergedunit)
write.csv(poecrAmbient_abund_mergedunit, "poecrAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
poecrAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(poecrAmbient_abund_mergedunit), data.frame(poecrAmbient_abund_mergedunit, row.names=NULL))
View(poecrAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poecrAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
poecrAmbient_abund_mergedunit_asdf<-as_data_frame(poecrAmbient_abund_mergedunit) 
View(poecrAmbient_abund_mergedunit_asdf)
#Sum rows
poecrAmbient_abund_mergedunit_asdf_total<-cbind(poecrAmbient_abund_mergedunit_asdf, total =rowSums(poecrAmbient_abund_mergedunit_asdf)) 
View(poecrAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "poecrAmbient_abund_mergedunit_rn_to_col" to "poecrAmbient_abund_mergedunit_rn_to_col_total"
poecrAmbient_abund_mergedunit_asdf_total$Unit<-poecrAmbient_abund_mergedunit_rn_to_col$Unit
View(poecrAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
poecrAmbient_abund_mergedunit_asdf_total_clean<-poecrAmbient_abund_mergedunit_asdf_total[,-(1:12)]
View(poecrAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
poecrAmbient_abund_mergedunit_asdf_total_clean<-cbind(poecrAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(poecrAmbient_abund_mergedunit_asdf_total_clean)
write.csv(poecrAmbient_abund_mergedunit_asdf_total_clean, "poecrAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poecrAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(poecrAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poecrAmbient_abund_mergedunit_asdf_total_clean_sum)

#Order==Tethyida
Ambientphyseqtethrmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Order=="Tethyida")
Richness_Ambientphyseqtethrmergedunit<-estimate_richness(Ambientphyseqtethrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqtethrmergedunit, "Richness_Ambientphyseqtethrmergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqtethrmergedunit)=sample_names(Ambientphyseqtethrmergedunit)
Richness_Ambientphyseqtethrmergedunit<-dfRowName(Richness_Ambientphyseqtethrmergedunit, name ="UNIT")
View(Richness_Ambientphyseqtethrmergedunit)
Richness_Ambientphyseqtethrmergedunit<-cbind(Richness_Ambientphyseqtethrmergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqtethrmergedunit_sum<-summarySE(Richness_Ambientphyseqtethrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqtethrmergedunit_sum)
#Shannon
Shannon_Ambientphyseqtethrmergedunit_sum<-summarySE(Richness_Ambientphyseqtethrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqtethrmergedunit_sum)
#Simpson
Simpson_Ambientphyseqtethrmergedunit_sum<-summarySE(Richness_Ambientphyseqtethrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqtethrmergedunit_sum)
#Getting abundance
tethrAmbient_abund_mergedunit=as(otu_table(Ambientphyseqtethrmergedunit), "matrix")
View(tethrAmbient_abund_mergedunit)
write.csv(tethrAmbient_abund_mergedunit, "tethrAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
tethrAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(tethrAmbient_abund_mergedunit), data.frame(tethrAmbient_abund_mergedunit, row.names=NULL))
View(tethrAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tethrAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
tethrAmbient_abund_mergedunit_asdf<-as_data_frame(tethrAmbient_abund_mergedunit) 
View(tethrAmbient_abund_mergedunit_asdf)
#Sum rows
tethrAmbient_abund_mergedunit_asdf_total<-cbind(tethrAmbient_abund_mergedunit_asdf, total =rowSums(tethrAmbient_abund_mergedunit_asdf)) 
View(tethrAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "tethrAmbient_abund_mergedunit_rn_to_col" to "tethrAmbient_abund_mergedunit_rn_to_col_total"
tethrAmbient_abund_mergedunit_asdf_total$Unit<-tethrAmbient_abund_mergedunit_rn_to_col$Unit
View(tethrAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
tethrAmbient_abund_mergedunit_asdf_total_clean<-tethrAmbient_abund_mergedunit_asdf_total[,-(1:6)]
View(tethrAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
tethrAmbient_abund_mergedunit_asdf_total_clean<-cbind(tethrAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(tethrAmbient_abund_mergedunit_asdf_total_clean)
write.csv(tethrAmbient_abund_mergedunit_asdf_total_clean, "tethrAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tethrAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(tethrAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tethrAmbient_abund_mergedunit_asdf_total_clean_sum)

#Order==Haplosclerida
Ambientphyseqhaprmergedunit<-subset_taxa(Ambientphyseqpormergedunit, Order=="Haplosclerida")
Richness_Ambientphyseqhaprmergedunit<-estimate_richness(Ambientphyseqhaprmergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_Ambientphyseqhaprmergedunit, "Richness_Ambientphyseqhaprmergedunit.csv", row.names = F)
row.names(Richness_Ambientphyseqhaprmergedunit)=sample_names(Ambientphyseqhaprmergedunit)
Richness_Ambientphyseqhaprmergedunit<-dfRowName(Richness_Ambientphyseqhaprmergedunit, name ="UNIT")
View(Richness_Ambientphyseqhaprmergedunit)
Richness_Ambientphyseqhaprmergedunit<-cbind(Richness_Ambientphyseqhaprmergedunit, TMT='ATACO2')
#Observed
Richness_Ambientphyseqhaprmergedunit_sum<-summarySE(Richness_Ambientphyseqhaprmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_Ambientphyseqhaprmergedunit_sum)
#Shannon
Shannon_Ambientphyseqhaprmergedunit_sum<-summarySE(Richness_Ambientphyseqhaprmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_Ambientphyseqhaprmergedunit_sum)
#Simpson
Simpson_Ambientphyseqhaprmergedunit_sum<-summarySE(Richness_Ambientphyseqhaprmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_Ambientphyseqhaprmergedunit_sum)
#Getting abundance
haprAmbient_abund_mergedunit=as(otu_table(Ambientphyseqhaprmergedunit), "matrix")
View(haprAmbient_abund_mergedunit)
write.csv(haprAmbient_abund_mergedunit, "haprAmbient_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
haprAmbient_abund_mergedunit_rn_to_col <- cbind(rownames(haprAmbient_abund_mergedunit), data.frame(haprAmbient_abund_mergedunit, row.names=NULL))
View(haprAmbient_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(haprAmbient_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
haprAmbient_abund_mergedunit_asdf<-as_data_frame(haprAmbient_abund_mergedunit) 
View(haprAmbient_abund_mergedunit_asdf)
#Sum rows
haprAmbient_abund_mergedunit_asdf_total<-cbind(haprAmbient_abund_mergedunit_asdf, total =rowSums(haprAmbient_abund_mergedunit_asdf)) 
View(haprAmbient_abund_mergedunit_asdf_total)
#add column "Unit" from "haprAmbient_abund_mergedunit_rn_to_col" to "haprAmbient_abund_mergedunit_rn_to_col_total"
haprAmbient_abund_mergedunit_asdf_total$Unit<-haprAmbient_abund_mergedunit_rn_to_col$Unit
names(haprAmbient_abund_mergedunit_asdf_total)
View(haprAmbient_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
haprAmbient_abund_mergedunit_asdf_total_clean<-haprAmbient_abund_mergedunit_asdf_total[,-(1:22)]
View(haprAmbient_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "Ambient"
haprAmbient_abund_mergedunit_asdf_total_clean<-cbind(haprAmbient_abund_mergedunit_asdf_total_clean, TMT='ATACO2')
View(haprAmbient_abund_mergedunit_asdf_total_clean)
write.csv(haprAmbient_abund_mergedunit_asdf_total_clean, "haprAmbient_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
haprAmbient_abund_mergedunit_asdf_total_clean_sum<-summarySE(haprAmbient_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(haprAmbient_abund_mergedunit_asdf_total_clean_sum)


#Rarefraction curve preliminary steps:Subset data within "ATACO2DATE"to obtain the different subsponge groups.
#remember that when subsetting taxonomic assignments you must use "subset_taxa" starting with all porifera

#Subset data by TMT="Ambient" 
Ambientphyseqpor=subset_samples(physeqpor,TMT=="ATACO2")
AmbientphyseqpormergedDate = merge_samples(Ambientphyseqpor, "Date")
SD= merge_samples(sample_data(Ambientphyseqpor), "Date")
print(SD[, "Date"])
print(AmbientphyseqpormergedDate)
sample_names(AmbientphyseqpormergedDate)
sample_data(AmbientphyseqpormergedDate)$Date
identical(SD,sample_data(AmbientphyseqpormergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqpormergedDate)$Date <- factor(sample_names(AmbientphyseqpormergedDate))
View(AmbientphyseqpormergedDate)
#Getting abundance
poriferaAmbient_abund_mergedDate=as(otu_table(AmbientphyseqpormergedDate), "matrix")
View(poriferaAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poriferaAmbient_abund_mergedDate_abs_pre<-poriferaAmbient_abund_mergedDate
write.csv(poriferaAmbient_abund_mergedDate, "poriferaAmbient_abund_mergedDate.csv", row.names = F)
poriferaAmbient_abund_mergedDate_abs_pre[poriferaAmbient_abund_mergedDate_abs_pre >0] <-1 #making a presence absence data frame
View(poriferaAmbient_abund_mergedDate_abs_pre)
write.csv(poriferaAmbient_abund_mergedDate_abs_pre, "poriferaAmbient_abund_mergedDate_abs_pre.csv", row.names = F)
#in excel I merged both ambient and intake .csv files and reimported into "R"

#Subset data within "Ambient" to obtain the different subsponge groups.
#remember that when subsetting taxonomic assignments you must use "subset_taxa" starting with 
#Homoscleromorpha
Ambientphyseqhomo=subset_taxa(Ambientphyseqpor,Class=="Homoscleromorpha")
AmbientphyseqhomomergedDate = merge_samples(Ambientphyseqhomo, "Date")
SD= merge_samples(sample_data(Ambientphyseqhomo), "Date")
print(SD[, "Date"])
print(AmbientphyseqhomomergedDate)
sample_names(AmbientphyseqhomomergedDate)
sample_data(AmbientphyseqhomomergedDate)$Date
identical(SD,sample_data(AmbientphyseqhomomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqhomomergedDate)$Date <- factor(sample_names(AmbientphyseqhomomergedDate))
#Getting abundance
homoAmbient_abund_mergedDate=as(otu_table(AmbientphyseqhomomergedDate), "matrix")
View(homoAmbient_abund_mergedDate)
#adding column at the end that shows TMT=Homoscleromorpha
homoAmbient_abund_mergedDate<-cbind(homoAmbient_abund_mergedDate, TMT='ATACO2')
head(homoAmbient_abund_mergedDate)
write.csv(homoAmbient_abund_mergedDate, "homoAmbient_abund_mergedDate.csv", row.names = F)

#Demospongiae
Ambientphyseqdem=subset_taxa(Ambientphyseqpor,Class=="Demospongiae")
AmbientphyseqdemmergedDate = merge_samples(Ambientphyseqdem, "Date")
SD= merge_samples(sample_data(Ambientphyseqdem), "Date")
print(SD[, "Date"])
print(AmbientphyseqdemmergedDate)
sample_names(AmbientphyseqdemmergedDate)
sample_data(AmbientphyseqdemmergedDate)$Date
identical(SD,sample_data(AmbientphyseqdemmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqdemmergedDate)$Date <- factor(sample_names(AmbientphyseqdemmergedDate))
#Getting abundance
demAmbient_abund_mergedDate=as(otu_table(AmbientphyseqdemmergedDate), "matrix")
View(demAmbient_abund_mergedDate)
#adding column at the end that shows TMT=ATACO2
demAmbient_abund_mergedDate<-cbind(demAmbient_abund_mergedDate, TMT='ATACO2')
write.csv(demAmbient_abund_mergedDate, "demAmbient_abund_mergedDate.csv", row.names = F)


###Subset with Clathrinida 
Ambientphyseqclat=subset_taxa(Ambientphyseqpor,Order=="Clathrinida")
AmbientphyseqclatrmergedDate = merge_samples(Ambientphyseqclat, "Date")
SD= merge_samples(sample_data(Ambientphyseqclat), "Date")
print(SD[, "Date"])
print(AmbientphyseqclatrmergedDate)
sample_names(AmbientphyseqclatrmergedDate)
sample_data(AmbientphyseqclatrmergedDate)$Date
identical(SD,sample_data(AmbientphyseqclatrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqclatrmergedDate)$Date <- factor(sample_names(AmbientphyseqclatrmergedDate))
#Getting abundance
clatAmbient_abund_mergedDate=as(otu_table(AmbientphyseqclatrmergedDate), "matrix")
View(clatAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
clatAmbient_abund_mergedDate_abs_pre<-clatAmbient_abund_mergedDate
write.csv(clatAmbient_abund_mergedDate, "clatAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Leucosolenida
Ambientphyseqleuc=subset_taxa(Ambientphyseqpor,Order=="Leucosolenida")
AmbientphyseqleucmergedDate = merge_samples(Ambientphyseqleuc, "Date")
SD= merge_samples(sample_data(Ambientphyseqleuc), "Date")
print(SD[, "Date"])
print(AmbientphyseqleucmergedDate)
sample_names(AmbientphyseqleucmergedDate)
sample_data(AmbientphyseqleucmergedDate)$Date
identical(SD,sample_data(AmbientphyseqleucmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqleucmergedDate)$Date <- factor(sample_names(AmbientphyseqleucmergedDate))
#Getting abundance
leucAmbient_abund_mergedDate=as(otu_table(AmbientphyseqleucmergedDate), "matrix")
View(leucAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
leucAmbient_abund_mergedDate_abs_pre<-leucAmbient_abund_mergedDate
write.csv(leucAmbient_abund_mergedDate, "leucAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Class: Calcarea
Ambientphyseqcalc=subset_taxa(Ambientphyseqpor,Class=="Calcarea")
AmbientphyseqcalcmergedDate = merge_samples(Ambientphyseqcalc, "Date")
SD= merge_samples(sample_data(Ambientphyseqcalc), "Date")
print(SD[, "Date"])
print(AmbientphyseqcalcmergedDate)
sample_names(AmbientphyseqcalcmergedDate)
sample_data(AmbientphyseqcalcmergedDate)$Date
identical(SD,sample_data(AmbientphyseqcalcmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqcalcmergedDate)$Date <- factor(sample_names(AmbientphyseqcalcmergedDate))
#Getting abundance
calcAmbient_abund_mergedDate=as(otu_table(AmbientphyseqcalcmergedDate), "matrix")
View(calcAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
calcAmbient_abund_mergedDate_abs_pre<-calcAmbient_abund_mergedDate
#adding column at the end that shows TMT=ATACO2
calcAmbient_abund_mergedDate<-cbind(calcAmbient_abund_mergedDate, TMT='ATACO2')
write.csv(calcAmbient_abund_mergedDate, "calcAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Keratosa 
Ambientphyseqker=subset_taxa(Ambientphyseqpor,Subclass=="Keratosa")
AmbientphyseqkermergedDate = merge_samples(Ambientphyseqker, "Date")
SD= merge_samples(sample_data(Ambientphyseqker), "Date")
print(SD[, "Date"])
print(AmbientphyseqkermergedDate)
sample_names(AmbientphyseqkermergedDate)
sample_data(AmbientphyseqkermergedDate)$Date
identical(SD,sample_data(AmbientphyseqkermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqkermergedDate)$Date <- factor(sample_names(AmbientphyseqkermergedDate))
#Getting abundance
kerAmbient_abund_mergedDate=as(otu_table(AmbientphyseqkermergedDate), "matrix")
View(kerAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
kerAmbient_abund_mergedDate_abs_pre<-kerAmbient_abund_mergedDate
write.csv(kerAmbient_abund_mergedDate, "kerAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Tetractinellida 
Ambientphyseqtetr=subset_taxa(Ambientphyseqpor,Order=="Tetractinellida")
AmbientphyseqtetrmergedDate = merge_samples(Ambientphyseqtetr, "Date")
SD= merge_samples(sample_data(Ambientphyseqtetr), "Date")
print(SD[, "Date"])
print(AmbientphyseqtetrmergedDate)
sample_names(AmbientphyseqtetrmergedDate)
sample_data(AmbientphyseqtetrmergedDate)$Date
identical(SD,sample_data(AmbientphyseqtetrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqtetrmergedDate)$Date <- factor(sample_names(AmbientphyseqtetrmergedDate))
#Getting abundance
tetrAmbient_abund_mergedDate=as(otu_table(AmbientphyseqtetrmergedDate), "matrix")
View(tetrAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tetrAmbient_abund_mergedDate_abs_pre<-tetrAmbient_abund_mergedDate
write.csv(tetrAmbient_abund_mergedDate, "tetrAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Suberitida 
Ambientphyseqsub=subset_taxa(Ambientphyseqpor,Order=="Suberitida")
AmbientphyseqsubmergedDate = merge_samples(Ambientphyseqsub, "Date")
SD= merge_samples(sample_data(Ambientphyseqsub), "Date")
print(SD[, "Date"])
print(AmbientphyseqsubmergedDate)
sample_names(AmbientphyseqsubmergedDate)
sample_data(AmbientphyseqsubmergedDate)$Date
identical(SD,sample_data(AmbientphyseqsubmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqsubmergedDate)$Date <- factor(sample_names(AmbientphyseqsubmergedDate))
#Getting abundance
subAmbient_abund_mergedDate=as(otu_table(AmbientphyseqsubmergedDate), "matrix")
View(subAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
subAmbient_abund_mergedDate_abs_pre<-subAmbient_abund_mergedDate
write.csv(subAmbient_abund_mergedDate, "subAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Poecilosclerida
Ambientphyseqpoec=subset_taxa(Ambientphyseqpor,Order=="Poecilosclerida")
AmbientphyseqpoecmergedDate = merge_samples(Ambientphyseqpoec, "Date")
SD= merge_samples(sample_data(Ambientphyseqpoec), "Date")
print(SD[, "Date"])
print(AmbientphyseqpoecmergedDate)
sample_names(AmbientphyseqpoecmergedDate)
sample_data(AmbientphyseqpoecmergedDate)$Date
identical(SD,sample_data(AmbientphyseqpoecmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqpoecmergedDate)$Date <- factor(sample_names(AmbientphyseqpoecmergedDate))
#Getting abundance
poecAmbient_abund_mergedDate=as(otu_table(AmbientphyseqpoecmergedDate), "matrix")
View(poecAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poecAmbient_abund_mergedDate_abs_pre<-poecAmbient_abund_mergedDate
write.csv(poecAmbient_abund_mergedDate, "poecAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Haplosclerida
Ambientphyseqhap=subset_taxa(Ambientphyseqpor,Order=="Haplosclerida")
AmbientphyseqhapmergedDate = merge_samples(Ambientphyseqhap, "Date")
SD= merge_samples(sample_data(Ambientphyseqhap), "Date")
print(SD[, "Date"])
print(AmbientphyseqhapmergedDate)
sample_names(AmbientphyseqhapmergedDate)
sample_data(AmbientphyseqhapmergedDate)$Date
identical(SD,sample_data(AmbientphyseqhapmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqhapmergedDate)$Date <- factor(sample_names(AmbientphyseqhapmergedDate))
#Getting abundance
hapAmbient_abund_mergedDate=as(otu_table(AmbientphyseqhapmergedDate), "matrix")
View(hapAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
hapAmbient_abund_mergedDate_abs_pre<-hapAmbient_abund_mergedDate
write.csv(hapAmbient_abund_mergedDate, "hapAmbient_abund_mergedDate.csv", row.names = F)

###Subset with Tethyida
Ambientphyseqteth=subset_taxa(Ambientphyseqpor,Order=="Tethyida")
AmbientphyseqtethmergedDate = merge_samples(Ambientphyseqteth, "Date")
SD= merge_samples(sample_data(Ambientphyseqteth), "Date")
print(SD[, "Date"])
print(AmbientphyseqtethmergedDate)
sample_names(AmbientphyseqtethmergedDate)
sample_data(AmbientphyseqtethmergedDate)$Date
identical(SD,sample_data(AmbientphyseqtethmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(AmbientphyseqtethmergedDate)$Date <- factor(sample_names(AmbientphyseqtethmergedDate))
#Getting abundance
tethAmbient_abund_mergedDate=as(otu_table(AmbientphyseqtethmergedDate), "matrix")
View(tethAmbient_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tethAmbient_abund_mergedDate_abs_pre<-tethAmbient_abund_mergedDate
write.csv(tethAmbient_abund_mergedDate, "tethAmbient_abund_mergedDate.csv", row.names = F)

#merge all of the csv per sponge group and bring in to R
poriferabyclasorder_ambient_abund_mergedDate<-read.csv("poriferabyclasorder_ambient_abund_mergedDate.csv") #import .csv
View(poriferabyclasorder_ambient_abund_mergedDate)
poriferabyclasorder_ambient_abund_mergedDate[is.na(poriferabyclasorder_ambient_abund_mergedDate)] <- 0
View(poriferabyclasorder_ambient_abund_mergedDate)
Porif_ambient_class_order<-poriferabyclasorder_ambient_abund_mergedDate
View(Porif_ambient_class_order)
names(Porif_ambient_class_order)
####rarefraction curve###
Porif_ambient_class_order_all<-Porif_ambient_class_order[2:94]
curve_Porif_ambient_class_order_all = specaccum(Porif_ambient_class_order_all, method = "rarefaction", 
                                                permutations = 100)
#subset each habitat into its own df
Porif_ambient_class_order%>% filter(Site == "Calcarea") -> Calcarea
names(Calcarea)
Porif_ambient_class_order%>% filter(Site == "Homoscleromorpha") -> Homoscleromorpha
Porif_ambient_class_order%>% filter(Site == "Keratosa") -> Keratosa
Porif_ambient_class_order%>% filter(Site == "Tetractinellida") -> Tetractinellida
Porif_ambient_class_order%>% filter(Site == "Suberitida") -> Suberitida
Porif_ambient_class_order%>% filter(Site == "Poecilosclerida") -> Poecilosclerida
Porif_ambient_class_order%>% filter(Site == "Haplosclerida") -> Haplosclerida
Porif_ambient_class_order%>% filter(Site == "Tethyida") -> Tethyida

#species accumulation curve for each habitat using all sponges
curve_amb_calc = specaccum(Calcarea[, 2:94], method = "rarefaction")
curve_amb_hom = specaccum(Homoscleromorpha[, 2:94], method = "rarefaction")
curve_amb_ker = specaccum(Keratosa[, 2:94], method = "rarefaction")
curve_amb_tetr = specaccum(Tetractinellida[, 2:94], method = "rarefaction")
curve_amb_sub = specaccum(Suberitida[, 2:94], method = "rarefaction")
curve_amb_poec = specaccum(Poecilosclerida[, 2:94], method = "rarefaction")
curve_amb_hap = specaccum(Haplosclerida[, 2:94], method = "rarefaction")
curve_amb_teth = specaccum(Tethyida[, 2:94], method = "rarefaction")

#figure out color palet "display.brewer.pal(n = 8, name = 'RdBu')"
#got the error Error in plot.new() : figure margins too large" so I typed the following command "par(mar=c(5,5,3,1))"" # this fixed the margins and shows axis vallues. 
#plot curve_all first by dates
par(mar=c(5,5,3,1))
plot(curve_amb_calc, xvar=c("individuals"),ylim=c(0,15),xlim=c(0,800),col="#B2182B", lwd=2, ci.lty=0, ci.col="#F4A582",
     main = "Default: Prettier CI")

#then plot the rest
plot(curve_amb_hom, add = TRUE,xvar=c("individuals"), col="#2166AC", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI") #col is COLOUR setting, so change it to something else if you 

plot(curve_amb_ker, add = TRUE,xvar=c("individuals"), col="#E1BE6A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_amb_tetr, add = TRUE,xvar=c("individuals"), col="#40B0A6", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_amb_sub, add = TRUE,xvar=c("individuals"), col="#E66100", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_amb_poec, add = TRUE,xvar=c("individuals"), col="#5D3A9B", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_amb_hap, add = TRUE,xvar=c("individuals"), col="#1AFF1A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_amb_teth, add = TRUE,xvar=c("individuals"), col="#994F00", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")
legend(500, 1, legend=c("Calcareous", "Homoscleromorpha"),
       col=c("#B2182B", "#2166AC"), lty=1:1, cex=0.8)
```

#Subset data by TMT="HTACO2" 
```{r setup, include=FALSE}
#Subset data by TMT="HTACO2" 
HTACO2physeqpor=subset_samples(physeqpor,TMT=="HTACO2")
#OTHER MS For ARMS vs ReEf paper phylogenetic analysis etc. 
HTACO2physeqpormergedDate = merge_samples(HTACO2physeqpor, "Date")
SD= merge_samples(sample_data(HTACO2physeqpor), "Date")
print(SD[, "Date"])
print(HTACO2physeqpormergedDate)
sample_names(HTACO2physeqpormergedDate)
sample_data(HTACO2physeqpormergedDate)$Date
identical(SD,sample_data(HTACO2physeqpormergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqpormergedDate)$Date <- factor(sample_names(HTACO2physeqpormergedDate))
#Getting abundance
poriferaHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqpormergedDate), "matrix")
View(poriferaHTACO2_abund_mergedDate)
write.csv(poriferaHTACO2_abund_mergedDate, "poriferaHTACO2_abund_mergedDate.csv", row.names = F)

#incorporating abundance changing it to presence absence data
poriferaHTACO2_abund_mergedDate_abs_pre<-poriferaHTACO2_abund_mergedDate

poriferaHTACO2_abund_mergedDate_abs_pre[poriferaHTACO2_abund_mergedDate_abs_pre >0] <-1 #making a presence absence data frame
View(poriferaHTACO2_abund_mergedDate_abs_pre)
write.csv(poriferaHTACO2_abund_mergedDate_abs_pre, "poriferaHTACO2_abund_mergedDate_abs_pre.csv", row.names = F)
#Getting abundance
poriferaHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqpormergedDate), "matrix")
View(poriferaHTACO2_abund_mergedunit)
write.csv(poriferaHTACO2_abund_mergedunit, "poriferaHTACO2_abund_mergedunit.csv", row.names = F)


###subsetting samples by "HTACO2"" and then merging sample by "Unit"" 
#which will combine the total diversity of the experiment for this treatment and calculate mean values of diversity#####
HTACO2physeqpormergedunit = merge_samples(HTACO2physeqpor, "UNIT")
SD= merge_samples(sample_data(HTACO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(HTACO2physeqpormergedunit)
sample_names(HTACO2physeqpormergedunit)
sample_data(HTACO2physeqpormergedunit)$UNIT
identical(SD,sample_data(HTACO2physeqpormergedunit))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqpormergedunit)$UNIT <- factor(sample_names(HTACO2physeqpormergedunit))

#Determine diversity within Porifera HTACO2.
Richness_HTACO2physeqpormergedunit<-estimate_richness(HTACO2physeqpormergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqpormergedunit, "Richness_HTACO2physeqpormergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqpormergedunit)=sample_names(HTACO2physeqpormergedunit)
Richness_HTACO2physeqpormergedunit<-dfRowName(Richness_HTACO2physeqpormergedunit, name ="UNIT")
View(Richness_HTACO2physeqpormergedunit)
#need to add a column "TMT" to indicate "HTACO2"
Richness_HTACO2physeqpormergedunit<-cbind(Richness_HTACO2physeqpormergedunit, TMT='HTACO2')
View(Richness_HTACO2physeqpormergedunit)
#Calculating mean values of Observed, SHannon, Simpson for the HTACO2. which will be shown in table xx diversity. 
#Observed
Richness_HTACO2physeqpormergedunit_sum<-summarySE(Richness_HTACO2physeqpormergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqpormergedunit_sum)
#Shannon
Shannon_HTACO2physeqpormergedunit_sum<-summarySE(Richness_HTACO2physeqpormergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqpormergedunit_sum)
#Simpson
Simpson_HTACO2physeqpormergedunit_sum<-summarySE(Richness_HTACO2physeqpormergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqpormergedunit_sum)

#Getting abundance
poriferaHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqpormergedunit), "matrix")
View(poriferaHTACO2_abund_mergedunit)
write.csv(poriferaHTACO2_abund_mergedunit, "poriferaHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
poriferaHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(poriferaHTACO2_abund_mergedunit), data.frame(poriferaHTACO2_abund_mergedunit, row.names=NULL))
View(poriferaHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poriferaHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
poriferaHTACO2_abund_mergedunit_asdf<-as_data_frame(poriferaHTACO2_abund_mergedunit) 
View(poriferaHTACO2_abund_mergedunit_asdf)
#Sum rows
poriferaHTACO2_abund_mergedunit_asdf_total<-cbind(poriferaHTACO2_abund_mergedunit_asdf, total =rowSums(poriferaHTACO2_abund_mergedunit_asdf)) 
View(poriferaHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poriferaHTACO2_abund_mergedunit_rn_to_col" to "poriferaHTACO2_abund_mergedunit_rn_to_col_total"
poriferaHTACO2_abund_mergedunit_asdf_total$Unit<-poriferaHTACO2_abund_mergedunit_rn_to_col$Unit
names(poriferaHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
poriferaHTACO2_abund_mergedunit_asdf_total_clean<-poriferaHTACO2_abund_mergedunit_asdf_total[,-(1:93)]
names(poriferaHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
poriferaHTACO2_abund_mergedunit_asdf_total_clean<-cbind(poriferaHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(poriferaHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(poriferaHTACO2_abund_mergedunit_asdf_total_clean, "poriferaHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poriferaHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poriferaHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poriferaHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Now calculate diversity within the different sponge subgroups within "HTACO2physeqpormergedunit"" starting with Class Homoscleromorpha

HTACO2physeqhomomergedunit<-subset_taxa(HTACO2physeqpormergedunit, Class=="Homoscleromorpha")
View(HTACO2physeqhomomergedunit)

Richness_HTACO2physeqhomomergedunit<-estimate_richness(HTACO2physeqhomomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTACO2physeqhomomergedunit)

#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqhomomergedunit, "Richness_HTACO2physeqhomomergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqhomomergedunit)=sample_names(HTACO2physeqhomomergedunit)
View(Richness_HTACO2physeqhomomergedunit)
Richness_HTACO2physeqhomomergedunit<-dfRowName(Richness_HTACO2physeqhomomergedunit, name ="UNIT")
View(Richness_HTACO2physeqhomomergedunit)
#Observed #add column with TMT name HTACO2
Richness_HTACO2physeqhomomergedunit <- Richness_HTACO2physeqhomomergedunit %>%
  add_column(TMT = "HTACO2")
Richness_HTACO2physeqhomomergedunit
Richness_HTACO2physeqhomomergedunit_sum<-summarySE(Richness_HTACO2physeqhomomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqhomomergedunit_sum)
#Shannon
Shannon_HTACO2physeqhomomergedunit_sum<-summarySE(Richness_HTACO2physeqhomomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqhomomergedunit_sum)
#Simpson
Simpson_HTACO2physeqhomomergedunit_sum<-summarySE(Richness_HTACO2physeqhomomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqhomomergedunit_sum)

#Getting abundance
homoHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqhomomergedunit), "matrix")
View(homoHTACO2_abund_mergedunit)
write.csv(homoHTACO2_abund_mergedunit, "homoHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
homoHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(homoHTACO2_abund_mergedunit), data.frame(homoHTACO2_abund_mergedunit, row.names=NULL))
View(homoHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(homoHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
homoHTACO2_abund_mergedunit_asdf<-as_data_frame(homoHTACO2_abund_mergedunit) 
View(homoHTACO2_abund_mergedunit_asdf)
#Sum rows
homoHTACO2_abund_mergedunit_asdf_total<-cbind(homoHTACO2_abund_mergedunit_asdf, total =rowSums(homoHTACO2_abund_mergedunit_asdf)) 
View(homoHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "homoHTACO2_abund_mergedunit_rn_to_col" to "homoHTACO2_abund_mergedunit_rn_to_col_total"
homoHTACO2_abund_mergedunit_asdf_total$Unit<-homoHTACO2_abund_mergedunit_rn_to_col$Unit
View(homoHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
homoHTACO2_abund_mergedunit_asdf_total_clean<-homoHTACO2_abund_mergedunit_asdf_total[,-(1:8)]
View(homoHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
homoHTACO2_abund_mergedunit_asdf_total_clean<-cbind(homoHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(homoHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(homoHTACO2_abund_mergedunit_asdf_total_clean, "homoHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
homoHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(homoHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(homoHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Demospongiae
#Now calculate diversity within the different sponge subgroups within "HTACO2physeqpormergedunit"" starting with Class Demospongiae

HTACO2physeqdemmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Class=="Demospongiae")
Richness_HTACO2physeqdemmergedunit<-estimate_richness(HTACO2physeqdemmergedunit, measures = c("Observed", "Shannon", "Simpson"))
Richness_HTACO2physeqdemmergedunit
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqdemmergedunit, "Richness_HTACO2physeqdemmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqdemmergedunit)=sample_names(HTACO2physeqdemmergedunit)
View(Richness_HTACO2physeqdemmergedunit)
Richness_HTACO2physeqdemmergedunit<-dfRowName(Richness_HTACO2physeqdemmergedunit, name ="UNIT")
View(Richness_HTACO2physeqdemmergedunit)
#Observed #add column with TMT name HTACO2
Richness_HTACO2physeqdemmergedunit <- Richness_HTACO2physeqdemmergedunit %>%
  add_column(TMT = "HTACO2")
Richness_HTACO2physeqdemmergedunit
Richness_HTACO2physeqdemmergedunit_sum<-summarySE(Richness_HTACO2physeqdemmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqdemmergedunit_sum)
#Shannon
Shannon_HTACO2physeqdemmergedunit_sum<-summarySE(Richness_HTACO2physeqdemmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqdemmergedunit_sum)
#Simpson
Simpson_HTACO2physeqdemmergedunit_sum<-summarySE(Richness_HTACO2physeqdemmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqdemmergedunit_sum)
#Getting abundance
demHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqdemmergedunit), "matrix")
View(demHTACO2_abund_mergedunit)
write.csv(demHTACO2_abund_mergedunit, "demHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
demHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(demHTACO2_abund_mergedunit), data.frame(demHTACO2_abund_mergedunit, row.names=NULL))
View(demHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(demHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
demHTACO2_abund_mergedunit_asdf<-as_data_frame(demHTACO2_abund_mergedunit) 

#Sum rows
demHTACO2_abund_mergedunit_asdf_total<-cbind(demHTACO2_abund_mergedunit_asdf, total =rowSums(demHTACO2_abund_mergedunit_asdf)) 
View(demHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "demHTACO2_abund_mergedunit_rn_to_col" to "demHTACO2_abund_mergedunit_rn_to_col_total"
demHTACO2_abund_mergedunit_asdf_total$Unit<-demHTACO2_abund_mergedunit_rn_to_col$Unit
names(demHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 66. 
demHTACO2_abund_mergedunit_asdf_total_clean<-demHTACO2_abund_mergedunit_asdf_total[,-(1:66)]
View(demHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
demHTACO2_abund_mergedunit_asdf_total_clean<-cbind(demHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(demHTACO2_abund_mergedunit_asdf_total_clean)
name_rows(demHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(demHTACO2_abund_mergedunit_asdf_total_clean, "demHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
demHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(demHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(demHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Dictyoceratida
#Now calculate diversity within the different sponge subgroups within "HTACO2physeqpormergedunit"" starting with Order Dictyoceratida
HTACO2physeqdictyomergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Dictyoceratida")
Richness_HTACO2physeqdictyomergedunit<-estimate_richness(HTACO2physeqdictyomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTACO2physeqdictyomergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqdictyomergedunit, "Richness_HTACO2physeqdictyomergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqdictyomergedunit)=sample_names(HTACO2physeqdictyomergedunit)
View(Richness_HTACO2physeqdictyomergedunit)
Richness_HTACO2physeqdictyomergedunit<-dfRowName(Richness_HTACO2physeqdictyomergedunit, name ="UNIT")
View(Richness_HTACO2physeqdictyomergedunit)
#Observed #add column with TMT name HTACO2
Richness_HTACO2physeqdictyomergedunit <- Richness_HTACO2physeqdictyomergedunit %>%
  add_column(TMT = "HTACO2")
Richness_HTACO2physeqdictyomergedunit
Richness_HTACO2physeqdictyomergedunit_sum<-summarySE(Richness_HTACO2physeqdictyomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqdictyomergedunit_sum)
#Shannon
Shannon_HTACO2physeqdictyomergedunit_sum<-summarySE(Richness_HTACO2physeqdictyomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqdictyomergedunit_sum)
#Simpson
Simpson_HTACO2physeqdictyomergedunit_sum<-summarySE(Richness_HTACO2physeqdictyomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqdictyomergedunit_sum)
#Getting abundance
dictyoHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqdictyomergedunit), "matrix")
View(dictyoHTACO2_abund_mergedunit)
write.csv(dictyoHTACO2_abund_mergedunit, "dictyoHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dictyoHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(dictyoHTACO2_abund_mergedunit), data.frame(dictyoHTACO2_abund_mergedunit, row.names=NULL))
View(dictyoHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dictyoHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dictyoHTACO2_abund_mergedunit_asdf<-as_data_frame(dictyoHTACO2_abund_mergedunit) 
View(dictyoHTACO2_abund_mergedunit_asdf)
#Sum rows
dictyoHTACO2_abund_mergedunit_asdf_total<-cbind(dictyoHTACO2_abund_mergedunit_asdf, total =rowSums(dictyoHTACO2_abund_mergedunit_asdf)) 
View(dictyoHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dictyoHTACO2_abund_mergedunit_rn_to_col" to "dictyoHTACO2_abund_mergedunit_rn_to_col_total"
dictyoHTACO2_abund_mergedunit_asdf_total$Unit<-dictyoHTACO2_abund_mergedunit_rn_to_col$Unit
View(dictyoHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dictyoHTACO2_abund_mergedunit_asdf_total_clean<-dictyoHTACO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dictyoHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
dictyoHTACO2_abund_mergedunit_asdf_total_clean<-cbind(dictyoHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(dictyoHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(dictyoHTACO2_abund_mergedunit_asdf_total_clean, "dictyoHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dictyoHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dictyoHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dictyoHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Dendroceratida
#Now calculate diversity within the different sponge subgroups within "HTACO2physeqpormergedunit"" starting with Order Dendroceratida
HTACO2physeqdendromergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Dendroceratida")
Richness_HTACO2physeqdendromergedunit<-estimate_richness(HTACO2physeqdendromergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTACO2physeqdendromergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqdendromergedunit, "Richness_HTACO2physeqdendromergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqdendromergedunit)=sample_names(HTACO2physeqdendromergedunit)
View(Richness_HTACO2physeqdendromergedunit)
Richness_HTACO2physeqdendromergedunit<-dfRowName(Richness_HTACO2physeqdendromergedunit, name ="UNIT")
View(Richness_HTACO2physeqdendromergedunit)
#Observed #add column with TMT name HTACO2
Richness_HTACO2physeqdendromergedunit <- Richness_HTACO2physeqdendromergedunit %>%
  add_column(TMT = "HTACO2")
Richness_HTACO2physeqdendromergedunit_sum<-summarySE(Richness_HTACO2physeqdendromergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqdendromergedunit_sum)
#Shannon
Shannon_HTACO2physeqdendromergedunit_sum<-summarySE(Richness_HTACO2physeqdendromergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqdendromergedunit_sum)
#Simpson
Simpson_HTACO2physeqdendromergedunit_sum<-summarySE(Richness_HTACO2physeqdendromergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqdendromergedunit_sum)
#Getting abundance
dendroHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqdendromergedunit), "matrix")
View(dendroHTACO2_abund_mergedunit)
write.csv(dendroHTACO2_abund_mergedunit, "dendroHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dendroHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(dendroHTACO2_abund_mergedunit), data.frame(dendroHTACO2_abund_mergedunit, row.names=NULL))
View(dendroHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dendroHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dendroHTACO2_abund_mergedunit_asdf<-as_data_frame(dendroHTACO2_abund_mergedunit) 
View(dendroHTACO2_abund_mergedunit_asdf)
#Sum rows
dendroHTACO2_abund_mergedunit_asdf_total<-cbind(dendroHTACO2_abund_mergedunit_asdf, total =rowSums(dendroHTACO2_abund_mergedunit_asdf)) 
View(dendroHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dendroHTACO2_abund_mergedunit_rn_to_col" to "dendroHTACO2_abund_mergedunit_rn_to_col_total"
dendroHTACO2_abund_mergedunit_asdf_total$Unit<-dendroHTACO2_abund_mergedunit_rn_to_col$Unit
View(dendroHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dendroHTACO2_abund_mergedunit_asdf_total_clean<-dendroHTACO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dendroHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
dendroHTACO2_abund_mergedunit_asdf_total_clean<-cbind(dendroHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(dendroHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(dendroHTACO2_abund_mergedunit_asdf_total_clean, "dendroHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dendroHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dendroHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dendroHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Calcarea
HTACO2physeqcalcmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Class=="Calcarea")
Richness_HTACO2physeqcalcmergedunit<-estimate_richness(HTACO2physeqcalcmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqcalcmergedunit, "Richness_HTACO2physeqcalcmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqcalcmergedunit)=sample_names(HTACO2physeqcalcmergedunit)
Richness_HTACO2physeqcalcmergedunit<-dfRowName(Richness_HTACO2physeqcalcmergedunit, name ="UNIT")
View(Richness_HTACO2physeqcalcmergedunit)
Richness_HTACO2physeqcalcmergedunit<-cbind(Richness_HTACO2physeqcalcmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqcalcmergedunit_sum<-summarySE(Richness_HTACO2physeqcalcmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqcalcmergedunit_sum)
#Shannon
Shannon_HTACO2physeqcalcmergedunit_sum<-summarySE(Richness_HTACO2physeqcalcmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqcalcmergedunit_sum)
#Simpson
Simpson_HTACO2physeqcalcmergedunit_sum<-summarySE(Richness_HTACO2physeqcalcmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqcalcmergedunit_sum)
#Getting abundance
calcHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqcalcmergedunit), "matrix")
View(calcHTACO2_abund_mergedunit)
write.csv(calcHTACO2_abund_mergedunit, "calcHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
calcHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(calcHTACO2_abund_mergedunit), data.frame(calcHTACO2_abund_mergedunit, row.names=NULL))
View(calcHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(calcHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
calcHTACO2_abund_mergedunit_asdf<-as_data_frame(calcHTACO2_abund_mergedunit) 
View(calcHTACO2_abund_mergedunit_asdf)
#Sum rows
calcHTACO2_abund_mergedunit_asdf_total<-cbind(calcHTACO2_abund_mergedunit_asdf, total =rowSums(calcHTACO2_abund_mergedunit_asdf)) 
View(calcHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "calcHTACO2_abund_mergedunit_rn_to_col" to "calcHTACO2_abund_mergedunit_rn_to_col_total"
calcHTACO2_abund_mergedunit_asdf_total$Unit<-calcHTACO2_abund_mergedunit_rn_to_col$Unit
View(calcHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
calcHTACO2_abund_mergedunit_asdf_total_clean<-calcHTACO2_abund_mergedunit_asdf_total[,-(1:19)]
View(calcHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
calcHTACO2_abund_mergedunit_asdf_total_clean<-cbind(calcHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(calcHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(calcHTACO2_abund_mergedunit_asdf_total_clean, "calcHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
calcHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(calcHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(calcHTACO2_abund_mergedunit_asdf_total_clean_sum)


#Keratosa
HTACO2physeqkermergedunit<-subset_taxa(HTACO2physeqpormergedunit, Subclass=="Keratosa")
Richness_HTACO2physeqkermergedunit<-estimate_richness(HTACO2physeqkermergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqkermergedunit, "Richness_HTACO2physeqkermergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqkermergedunit)=sample_names(HTACO2physeqkermergedunit)
Richness_HTACO2physeqkermergedunit<-dfRowName(Richness_HTACO2physeqkermergedunit, name ="UNIT")
View(Richness_HTACO2physeqkermergedunit)
Richness_HTACO2physeqkermergedunit<-cbind(Richness_HTACO2physeqkermergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqkermergedunit_sum<-summarySE(Richness_HTACO2physeqkermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqkermergedunit_sum)
#Shannon
Shannon_HTACO2physeqkermergedunit_sum<-summarySE(Richness_HTACO2physeqkermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqkermergedunit_sum)
#Simpson
Simpson_HTACO2physeqkermergedunit_sum<-summarySE(Richness_HTACO2physeqkermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqkermergedunit_sum)
#Getting abundance
kerHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqkermergedunit), "matrix")
View(kerHTACO2_abund_mergedunit)
write.csv(kerHTACO2_abund_mergedunit, "kerHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
kerHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(kerHTACO2_abund_mergedunit), data.frame(kerHTACO2_abund_mergedunit, row.names=NULL))
View(kerHTACO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(kerHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
kerHTACO2_abund_mergedunit_asdf<-as_data_frame(kerHTACO2_abund_mergedunit) 
View(kerHTACO2_abund_mergedunit_asdf)
#Sum rows
kerHTACO2_abund_mergedunit_asdf_total<-cbind(kerHTACO2_abund_mergedunit_asdf, total =rowSums(kerHTACO2_abund_mergedunit_asdf)) 
View(kerHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "kerHTACO2_abund_mergedunit_rn_to_col" to "kerHTACO2_abund_mergedunit_rn_to_col_total"
kerHTACO2_abund_mergedunit_asdf_total$Unit<-kerHTACO2_abund_mergedunit_rn_to_col$Unit
View(kerHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
kerHTACO2_abund_mergedunit_asdf_total_clean<-kerHTACO2_abund_mergedunit_asdf_total[,-(1:6)]
View(kerHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
kerHTACO2_abund_mergedunit_asdf_total_clean<-cbind(kerHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(kerHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(kerHTACO2_abund_mergedunit_asdf_total_clean, "kerHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
kerHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(kerHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(kerHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Subclass==Verongimorpha
HTACO2physeqvermergedunit<-subset_taxa(HTACO2physeqpormergedunit, Subclass=="Verongimorpha")
Richness_HTACO2physeqvermergedunit<-estimate_richness(HTACO2physeqvermergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_HTACO2physeqvermergedunit, "Richness_HTACO2physeqvermergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqvermergedunit)=sample_names(HTACO2physeqvermergedunit)
Richness_HTACO2physeqvermergedunit<-dfRowName(Richness_HTACO2physeqvermergedunit, name ="UNIT")
View(Richness_HTACO2physeqvermergedunit)
Richness_HTACO2physeqvermergedunit<-cbind(Richness_HTACO2physeqvermergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqvermergedunit_sum<-summarySE(Richness_HTACO2physeqvermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqvermergedunit_sum)
#Shannon
Shannon_HTACO2physeqvermergedunit_sum<-summarySE(Richness_HTACO2physeqvermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqvermergedunit_sum)
#Simpson
Simpson_HTACO2physeqvermergedunit_sum<-summarySE(Richness_HTACO2physeqvermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqvermergedunit_sum)
#Getting abundance
verHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqvermergedunit), "matrix")
View(verHTACO2_abund_mergedunit)
write.csv(verHTACO2_abund_mergedunit, "verHTACO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
verHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(verHTACO2_abund_mergedunit), data.frame(verHTACO2_abund_mergedunit, row.names=NULL))
View(verHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(verHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
verHTACO2_abund_mergedunit_asdf<-as_data_frame(verHTACO2_abund_mergedunit) 
View(verHTACO2_abund_mergedunit_asdf)
#Sum rows
verHTACO2_abund_mergedunit_asdf_total<-cbind(verHTACO2_abund_mergedunit_asdf, total =rowSums(verHTACO2_abund_mergedunit_asdf)) 
View(verHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "verHTACO2_abund_mergedunit_rn_to_col" to "verHTACO2_abund_mergedunit_rn_to_col_total"
verHTACO2_abund_mergedunit_asdf_total$Unit<-verHTACO2_abund_mergedunit_rn_to_col$Unit
View(verHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
verHTACO2_abund_mergedunit_asdf_total_clean<-verHTACO2_abund_mergedunit_asdf_total[,-(1:2)]
View(verHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
verHTACO2_abund_mergedunit_asdf_total_clean<-cbind(verHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(verHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(verHTACO2_abund_mergedunit_asdf_total_clean, "verHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
verHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(verHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(verHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Suberitida
HTACO2physeqsubmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Suberitida")
Richness_HTACO2physeqsubmergedunit<-estimate_richness(HTACO2physeqsubmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqsubmergedunit, "Richness_HTACO2physeqsubmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqsubmergedunit)=sample_names(HTACO2physeqsubmergedunit)
Richness_HTACO2physeqsubmergedunit<-dfRowName(Richness_HTACO2physeqsubmergedunit, name ="UNIT")
View(Richness_HTACO2physeqsubmergedunit)
Richness_HTACO2physeqsubmergedunit<-cbind(Richness_HTACO2physeqsubmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqsubmergedunit_sum<-summarySE(Richness_HTACO2physeqsubmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqsubmergedunit_sum)
#Shannon
Shannon_HTACO2physeqsubmergedunit_sum<-summarySE(Richness_HTACO2physeqsubmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqsubmergedunit_sum)
#Simpson
Simpson_HTACO2physeqsubmergedunit_sum<-summarySE(Richness_HTACO2physeqsubmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqsubmergedunit_sum)
#Getting abundance
subHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqsubmergedunit), "matrix")
View(subHTACO2_abund_mergedunit)
write.csv(subHTACO2_abund_mergedunit, "subHTACO2_abund_mergedunit.csv", row.names = F)
#Consubt Rowname into first column 
subHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(subHTACO2_abund_mergedunit), data.frame(subHTACO2_abund_mergedunit, row.names=NULL))
View(subHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(subHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#consubt to dataframe
subHTACO2_abund_mergedunit_asdf<-as_data_frame(subHTACO2_abund_mergedunit) 
View(subHTACO2_abund_mergedunit_asdf)
#Sum rows
subHTACO2_abund_mergedunit_asdf_total<-cbind(subHTACO2_abund_mergedunit_asdf, total =rowSums(subHTACO2_abund_mergedunit_asdf)) 
View(subHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "subHTACO2_abund_mergedunit_rn_to_col" to "subHTACO2_abund_mergedunit_rn_to_col_total"
subHTACO2_abund_mergedunit_asdf_total$Unit<-subHTACO2_abund_mergedunit_rn_to_col$Unit
View(subHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
subHTACO2_abund_mergedunit_asdf_total_clean<-subHTACO2_abund_mergedunit_asdf_total[,-(1:15)]
View(subHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
subHTACO2_abund_mergedunit_asdf_total_clean<-cbind(subHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(subHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(subHTACO2_abund_mergedunit_asdf_total_clean, "subHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
subHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(subHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(subHTACO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Tetractinellida
HTACO2physeqtetrmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Tetractinellida")
Richness_HTACO2physeqtetrmergedunit<-estimate_richness(HTACO2physeqtetrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqtetrmergedunit, "Richness_HTACO2physeqtetrmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqtetrmergedunit)=sample_names(HTACO2physeqtetrmergedunit)
Richness_HTACO2physeqtetrmergedunit<-dfRowName(Richness_HTACO2physeqtetrmergedunit, name ="UNIT")
View(Richness_HTACO2physeqtetrmergedunit)
Richness_HTACO2physeqtetrmergedunit<-cbind(Richness_HTACO2physeqtetrmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqtetrmergedunit_sum<-summarySE(Richness_HTACO2physeqtetrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqtetrmergedunit_sum)
#Shannon
Shannon_HTACO2physeqtetrmergedunit_sum<-summarySE(Richness_HTACO2physeqtetrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqtetrmergedunit_sum)
#Simpson
Simpson_HTACO2physeqtetrmergedunit_sum<-summarySE(Richness_HTACO2physeqtetrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqtetrmergedunit_sum)
#Getting abundance
tetrHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqtetrmergedunit), "matrix")
View(tetrHTACO2_abund_mergedunit)
write.csv(tetrHTACO2_abund_mergedunit, "tetrHTACO2_abund_mergedunit.csv", row.names = F)
#Contetrt Rowname into first column 
tetrHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(tetrHTACO2_abund_mergedunit), data.frame(tetrHTACO2_abund_mergedunit, row.names=NULL))
View(tetrHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tetrHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contetrt to dataframe
tetrHTACO2_abund_mergedunit_asdf<-as_data_frame(tetrHTACO2_abund_mergedunit) 
View(tetrHTACO2_abund_mergedunit_asdf)
#Sum rows
tetrHTACO2_abund_mergedunit_asdf_total<-cbind(tetrHTACO2_abund_mergedunit_asdf, total =rowSums(tetrHTACO2_abund_mergedunit_asdf)) 
View(tetrHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tetrHTACO2_abund_mergedunit_rn_to_col" to "tetrHTACO2_abund_mergedunit_rn_to_col_total"
tetrHTACO2_abund_mergedunit_asdf_total$Unit<-tetrHTACO2_abund_mergedunit_rn_to_col$Unit
View(tetrHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 3 
tetrHTACO2_abund_mergedunit_asdf_total_clean<-tetrHTACO2_abund_mergedunit_asdf_total[,-(1:3)]
View(tetrHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
tetrHTACO2_abund_mergedunit_asdf_total_clean<-cbind(tetrHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(tetrHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(tetrHTACO2_abund_mergedunit_asdf_total_clean, "tetrHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tetrHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tetrHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tetrHTACO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Poecilosclerida
HTACO2physeqpoecrmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Poecilosclerida")
Richness_HTACO2physeqpoecrmergedunit<-estimate_richness(HTACO2physeqpoecrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqpoecrmergedunit, "Richness_HTACO2physeqpoecrmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqpoecrmergedunit)=sample_names(HTACO2physeqpoecrmergedunit)
Richness_HTACO2physeqpoecrmergedunit<-dfRowName(Richness_HTACO2physeqpoecrmergedunit, name ="UNIT")
View(Richness_HTACO2physeqpoecrmergedunit)
Richness_HTACO2physeqpoecrmergedunit<-cbind(Richness_HTACO2physeqpoecrmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTACO2physeqpoecrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqpoecrmergedunit_sum)
#Shannon
Shannon_HTACO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTACO2physeqpoecrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqpoecrmergedunit_sum)
#Simpson
Simpson_HTACO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTACO2physeqpoecrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqpoecrmergedunit_sum)
#Getting abundance
poecrHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqpoecrmergedunit), "matrix")
View(poecrHTACO2_abund_mergedunit)
write.csv(poecrHTACO2_abund_mergedunit, "poecrHTACO2_abund_mergedunit.csv", row.names = F)
#Conpoecrt Rowname into first column 
poecrHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(poecrHTACO2_abund_mergedunit), data.frame(poecrHTACO2_abund_mergedunit, row.names=NULL))
View(poecrHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poecrHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conpoecrt to dataframe
poecrHTACO2_abund_mergedunit_asdf<-as_data_frame(poecrHTACO2_abund_mergedunit) 
View(poecrHTACO2_abund_mergedunit_asdf)
#Sum rows
poecrHTACO2_abund_mergedunit_asdf_total<-cbind(poecrHTACO2_abund_mergedunit_asdf, total =rowSums(poecrHTACO2_abund_mergedunit_asdf)) 
View(poecrHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poecrHTACO2_abund_mergedunit_rn_to_col" to "poecrHTACO2_abund_mergedunit_rn_to_col_total"
poecrHTACO2_abund_mergedunit_asdf_total$Unit<-poecrHTACO2_abund_mergedunit_rn_to_col$Unit
View(poecrHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
poecrHTACO2_abund_mergedunit_asdf_total_clean<-poecrHTACO2_abund_mergedunit_asdf_total[,-(1:12)]
View(poecrHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
poecrHTACO2_abund_mergedunit_asdf_total_clean<-cbind(poecrHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(poecrHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(poecrHTACO2_abund_mergedunit_asdf_total_clean, "poecrHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poecrHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poecrHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poecrHTACO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Tethyida
HTACO2physeqtethrmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Tethyida")
Richness_HTACO2physeqtethrmergedunit<-estimate_richness(HTACO2physeqtethrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqtethrmergedunit, "Richness_HTACO2physeqtethrmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqtethrmergedunit)=sample_names(HTACO2physeqtethrmergedunit)
Richness_HTACO2physeqtethrmergedunit<-dfRowName(Richness_HTACO2physeqtethrmergedunit, name ="UNIT")
View(Richness_HTACO2physeqtethrmergedunit)
Richness_HTACO2physeqtethrmergedunit<-cbind(Richness_HTACO2physeqtethrmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqtethrmergedunit_sum<-summarySE(Richness_HTACO2physeqtethrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqtethrmergedunit_sum)
#Shannon
Shannon_HTACO2physeqtethrmergedunit_sum<-summarySE(Richness_HTACO2physeqtethrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqtethrmergedunit_sum)
#Simpson
Simpson_HTACO2physeqtethrmergedunit_sum<-summarySE(Richness_HTACO2physeqtethrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqtethrmergedunit_sum)
#Getting abundance
tethrHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqtethrmergedunit), "matrix")
View(tethrHTACO2_abund_mergedunit)
write.csv(tethrHTACO2_abund_mergedunit, "tethrHTACO2_abund_mergedunit.csv", row.names = F)
#Contethrt Rowname into first column 
tethrHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(tethrHTACO2_abund_mergedunit), data.frame(tethrHTACO2_abund_mergedunit, row.names=NULL))
View(tethrHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tethrHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contethrt to dataframe
tethrHTACO2_abund_mergedunit_asdf<-as_data_frame(tethrHTACO2_abund_mergedunit) 
View(tethrHTACO2_abund_mergedunit_asdf)
#Sum rows
tethrHTACO2_abund_mergedunit_asdf_total<-cbind(tethrHTACO2_abund_mergedunit_asdf, total =rowSums(tethrHTACO2_abund_mergedunit_asdf)) 
View(tethrHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tethrHTACO2_abund_mergedunit_rn_to_col" to "tethrHTACO2_abund_mergedunit_rn_to_col_total"
tethrHTACO2_abund_mergedunit_asdf_total$Unit<-tethrHTACO2_abund_mergedunit_rn_to_col$Unit
View(tethrHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
tethrHTACO2_abund_mergedunit_asdf_total_clean<-tethrHTACO2_abund_mergedunit_asdf_total[,-(1:6)]
View(tethrHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
tethrHTACO2_abund_mergedunit_asdf_total_clean<-cbind(tethrHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(tethrHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(tethrHTACO2_abund_mergedunit_asdf_total_clean, "tethrHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tethrHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tethrHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tethrHTACO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Haplosclerida
HTACO2physeqhaprmergedunit<-subset_taxa(HTACO2physeqpormergedunit, Order=="Haplosclerida")
Richness_HTACO2physeqhaprmergedunit<-estimate_richness(HTACO2physeqhaprmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTACO2physeqhaprmergedunit, "Richness_HTACO2physeqhaprmergedunit.csv", row.names = F)
row.names(Richness_HTACO2physeqhaprmergedunit)=sample_names(HTACO2physeqhaprmergedunit)
Richness_HTACO2physeqhaprmergedunit<-dfRowName(Richness_HTACO2physeqhaprmergedunit, name ="UNIT")
View(Richness_HTACO2physeqhaprmergedunit)
Richness_HTACO2physeqhaprmergedunit<-cbind(Richness_HTACO2physeqhaprmergedunit, TMT='HTACO2')
#Observed
Richness_HTACO2physeqhaprmergedunit_sum<-summarySE(Richness_HTACO2physeqhaprmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTACO2physeqhaprmergedunit_sum)
#Shannon
Shannon_HTACO2physeqhaprmergedunit_sum<-summarySE(Richness_HTACO2physeqhaprmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTACO2physeqhaprmergedunit_sum)
#Simpson
Simpson_HTACO2physeqhaprmergedunit_sum<-summarySE(Richness_HTACO2physeqhaprmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTACO2physeqhaprmergedunit_sum)
#Getting abundance
haprHTACO2_abund_mergedunit=as(otu_table(HTACO2physeqhaprmergedunit), "matrix")
View(haprHTACO2_abund_mergedunit)
write.csv(haprHTACO2_abund_mergedunit, "haprHTACO2_abund_mergedunit.csv", row.names = F)
#Conhaprt Rowname into first column 
haprHTACO2_abund_mergedunit_rn_to_col <- cbind(rownames(haprHTACO2_abund_mergedunit), data.frame(haprHTACO2_abund_mergedunit, row.names=NULL))
View(haprHTACO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(haprHTACO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conhaprt to dataframe
haprHTACO2_abund_mergedunit_asdf<-as_data_frame(haprHTACO2_abund_mergedunit) 
View(haprHTACO2_abund_mergedunit_asdf)
#Sum rows
haprHTACO2_abund_mergedunit_asdf_total<-cbind(haprHTACO2_abund_mergedunit_asdf, total =rowSums(haprHTACO2_abund_mergedunit_asdf)) 
View(haprHTACO2_abund_mergedunit_asdf_total)
#add column "Unit" from "haprHTACO2_abund_mergedunit_rn_to_col" to "haprHTACO2_abund_mergedunit_rn_to_col_total"
haprHTACO2_abund_mergedunit_asdf_total$Unit<-haprHTACO2_abund_mergedunit_rn_to_col$Unit
names(haprHTACO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
haprHTACO2_abund_mergedunit_asdf_total_clean<-haprHTACO2_abund_mergedunit_asdf_total[,-(1:22)]
View(haprHTACO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTACO2"
haprHTACO2_abund_mergedunit_asdf_total_clean<-cbind(haprHTACO2_abund_mergedunit_asdf_total_clean, TMT='HTACO2')
View(haprHTACO2_abund_mergedunit_asdf_total_clean)
write.csv(haprHTACO2_abund_mergedunit_asdf_total_clean, "haprHTACO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
haprHTACO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(haprHTACO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(haprHTACO2_abund_mergedunit_asdf_total_clean_sum)


#Rarefraction curve preliminary steps:Subset data within "HTACO2DATE"to obtain the different subsponge groups.
#remember that when subsetting taxonomic assignments you must use "subset_taxa" starting with 
#Homoscleromorpha
HTACO2physeqhomo=subset_taxa(HTACO2physeqpor,Class=="Homoscleromorpha")
HTACO2physeqhomomergedDate = merge_samples(HTACO2physeqhomo, "Date")
SD= merge_samples(sample_data(HTACO2physeqhomo), "Date")
print(SD[, "Date"])
print(HTACO2physeqhomomergedDate)
sample_names(HTACO2physeqhomomergedDate)
sample_data(HTACO2physeqhomomergedDate)$Date
identical(SD,sample_data(HTACO2physeqhomomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqhomomergedDate)$Date <- factor(sample_names(HTACO2physeqhomomergedDate))
#Getting abundance
homoHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqhomomergedDate), "matrix")
View(homoHTACO2_abund_mergedDate)
#adding column at the end that shows TMT=HTACO2
homoHTACO2_abund_mergedDate<-cbind(homoHTACO2_abund_mergedDate, TMT='HTACO2')
head(homoHTACO2_abund_mergedDate)
write.csv(homoHTACO2_abund_mergedDate, "homoHTACO2_abund_mergedDate.csv", row.names = F)

#Demospongiae
HTACO2physeqdem=subset_taxa(HTACO2physeqpor,Class=="Demospongiae")
HTACO2physeqdemmergedDate = merge_samples(HTACO2physeqdem, "Date")
SD= merge_samples(sample_data(HTACO2physeqdem), "Date")
print(SD[, "Date"])
print(HTACO2physeqdemmergedDate)
sample_names(HTACO2physeqdemmergedDate)
sample_data(HTACO2physeqdemmergedDate)$Date
identical(SD,sample_data(HTACO2physeqdemmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqdemmergedDate)$Date <- factor(sample_names(HTACO2physeqdemmergedDate))
#Getting abundance
demHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqdemmergedDate), "matrix")
#adding column at the end that shows TMT=HTACO2
demHTACO2_abund_mergedDate<-cbind(demHTACO2_abund_mergedDate, TMT='HTACO2')
head(demHTACO2_abund_mergedDate)
names(demHTACO2_abund_mergedDate)
write.csv(demHTACO2_abund_mergedDate, "demHTACO2_abund_mergedDate.csv", row.names = F)


#Dictyoceratida
HTACO2physeqdictyo=subset_taxa(HTACO2physeqpor,Order=="Dictyoceratida")
HTACO2physeqdictyomergedDate = merge_samples(HTACO2physeqdictyo, "Date")
SD= merge_samples(sample_data(HTACO2physeqdictyo), "Date")
print(SD[, "Date"])
print(HTACO2physeqdictyomergedDate)
sample_names(HTACO2physeqdictyomergedDate)
sample_data(HTACO2physeqdictyomergedDate)$Date
identical(SD,sample_data(HTACO2physeqdictyomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqdictyomergedDate)$Date <- factor(sample_names(HTACO2physeqdictyomergedDate))
#Getting abundance
dictyoHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqdictyomergedDate), "matrix")
View(dictyoHTACO2_abund_mergedDate)
write.csv(dictyoHTACO2_abund_mergedDate, "dictyoHTACO2_abund_mergedDate.csv", row.names = F)

#Dendroceratida
HTACO2physeqdendro=subset_taxa(HTACO2physeqpor,Order=="Dendroceratida")
HTACO2physeqdendromergedDate = merge_samples(HTACO2physeqdendro, "Date")
SD= merge_samples(sample_data(HTACO2physeqdendro), "Date")
print(SD[, "Date"])
print(HTACO2physeqdendromergedDate)
sample_names(HTACO2physeqdendromergedDate)
sample_data(HTACO2physeqdendromergedDate)$Date
identical(SD,sample_data(HTACO2physeqdendromergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqdendromergedDate)$Date <- factor(sample_names(HTACO2physeqdendromergedDate))
#Getting abundance
dendroHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqdendromergedDate), "matrix")
View(dendroHTACO2_abund_mergedDate)
write.csv(dendroHTACO2_abund_mergedDate, "dendroHTACO2_abund_mergedDate.csv", row.names = F)


###Subset with Clathrinida 
HTACO2physeqclat=subset_taxa(HTACO2physeqpor,Order=="Clathrinida")
HTACO2physeqclatrmergedDate = merge_samples(HTACO2physeqclat, "Date")
SD= merge_samples(sample_data(HTACO2physeqclat), "Date")
print(SD[, "Date"])
print(HTACO2physeqclatrmergedDate)
sample_names(HTACO2physeqclatrmergedDate)
sample_data(HTACO2physeqclatrmergedDate)$Date
identical(SD,sample_data(HTACO2physeqclatrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqclatrmergedDate)$Date <- factor(sample_names(HTACO2physeqclatrmergedDate))
#Getting abundance
clatHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqclatrmergedDate), "matrix")
View(clatHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
clatHTACO2_abund_mergedDate_abs_pre<-clatHTACO2_abund_mergedDate
write.csv(clatHTACO2_abund_mergedDate, "clatHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Leucosolenida
HTACO2physeqleuc=subset_taxa(HTACO2physeqpor,Order=="Leucosolenida")
HTACO2physeqleucmergedDate = merge_samples(HTACO2physeqleuc, "Date")
SD= merge_samples(sample_data(HTACO2physeqleuc), "Date")
print(SD[, "Date"])
print(HTACO2physeqleucmergedDate)
sample_names(HTACO2physeqleucmergedDate)
sample_data(HTACO2physeqleucmergedDate)$Date
identical(SD,sample_data(HTACO2physeqleucmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqleucmergedDate)$Date <- factor(sample_names(HTACO2physeqleucmergedDate))
#Getting abundance
leucHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqleucmergedDate), "matrix")
View(leucHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
leucHTACO2_abund_mergedDate_abs_pre<-leucHTACO2_abund_mergedDate
write.csv(leucHTACO2_abund_mergedDate, "leucHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Class: Calcarea
HTACO2physeqcalc=subset_taxa(HTACO2physeqpor,Class=="Calcarea")
HTACO2physeqcalcmergedDate = merge_samples(HTACO2physeqcalc, "Date")
SD= merge_samples(sample_data(HTACO2physeqcalc), "Date")
print(SD[, "Date"])
print(HTACO2physeqcalcmergedDate)
sample_names(HTACO2physeqcalcmergedDate)
sample_data(HTACO2physeqcalcmergedDate)$Date
identical(SD,sample_data(HTACO2physeqcalcmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqcalcmergedDate)$Date <- factor(sample_names(HTACO2physeqcalcmergedDate))
#Getting abundance
calcHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqcalcmergedDate), "matrix")
View(calcHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
calcHTACO2_abund_mergedDate_abs_pre<-calcHTACO2_abund_mergedDate
#adding column at the end that shows TMT=HTACO2
calcHTACO2_abund_mergedDate<-cbind(calcHTACO2_abund_mergedDate, TMT='HTACO2')
write.csv(calcHTACO2_abund_mergedDate, "calcHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Keratosa 
HTACO2physeqker=subset_taxa(HTACO2physeqpor,Subclass=="Keratosa")
HTACO2physeqkermergedDate = merge_samples(HTACO2physeqker, "Date")
SD= merge_samples(sample_data(HTACO2physeqker), "Date")
print(SD[, "Date"])
print(HTACO2physeqkermergedDate)
sample_names(HTACO2physeqkermergedDate)
sample_data(HTACO2physeqkermergedDate)$Date
identical(SD,sample_data(HTACO2physeqkermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqkermergedDate)$Date <- factor(sample_names(HTACO2physeqkermergedDate))
#Getting abundance
kerHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqkermergedDate), "matrix")
View(kerHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
kerHTACO2_abund_mergedDate_abs_pre<-kerHTACO2_abund_mergedDate
write.csv(kerHTACO2_abund_mergedDate, "kerHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Tetractinellida 
HTACO2physeqtetr=subset_taxa(HTACO2physeqpor,Order=="Tetractinellida")
HTACO2physeqtetrmergedDate = merge_samples(HTACO2physeqtetr, "Date")
SD= merge_samples(sample_data(HTACO2physeqtetr), "Date")
print(SD[, "Date"])
print(HTACO2physeqtetrmergedDate)
sample_names(HTACO2physeqtetrmergedDate)
sample_data(HTACO2physeqtetrmergedDate)$Date
identical(SD,sample_data(HTACO2physeqtetrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqtetrmergedDate)$Date <- factor(sample_names(HTACO2physeqtetrmergedDate))
#Getting abundance
tetrHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqtetrmergedDate), "matrix")
View(tetrHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tetrHTACO2_abund_mergedDate_abs_pre<-tetrHTACO2_abund_mergedDate
write.csv(tetrHTACO2_abund_mergedDate, "tetrHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Suberitida 
HTACO2physeqsub=subset_taxa(HTACO2physeqpor,Order=="Suberitida")
HTACO2physeqsubmergedDate = merge_samples(HTACO2physeqsub, "Date")
SD= merge_samples(sample_data(HTACO2physeqsub), "Date")
print(SD[, "Date"])
print(HTACO2physeqsubmergedDate)
sample_names(HTACO2physeqsubmergedDate)
sample_data(HTACO2physeqsubmergedDate)$Date
identical(SD,sample_data(HTACO2physeqsubmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqsubmergedDate)$Date <- factor(sample_names(HTACO2physeqsubmergedDate))
#Getting abundance
subHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqsubmergedDate), "matrix")
View(subHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
subHTACO2_abund_mergedDate_abs_pre<-subHTACO2_abund_mergedDate
write.csv(subHTACO2_abund_mergedDate, "subHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Poecilosclerida
HTACO2physeqpoec=subset_taxa(HTACO2physeqpor,Order=="Poecilosclerida")
HTACO2physeqpoecmergedDate = merge_samples(HTACO2physeqpoec, "Date")
SD= merge_samples(sample_data(HTACO2physeqpoec), "Date")
print(SD[, "Date"])
print(HTACO2physeqpoecmergedDate)
sample_names(HTACO2physeqpoecmergedDate)
sample_data(HTACO2physeqpoecmergedDate)$Date
identical(SD,sample_data(HTACO2physeqpoecmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqpoecmergedDate)$Date <- factor(sample_names(HTACO2physeqpoecmergedDate))
#Getting abundance
poecHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqpoecmergedDate), "matrix")
View(poecHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poecHTACO2_abund_mergedDate_abs_pre<-poecHTACO2_abund_mergedDate
write.csv(poecHTACO2_abund_mergedDate, "poecHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Haplosclerida
HTACO2physeqhap=subset_taxa(HTACO2physeqpor,Order=="Haplosclerida")
HTACO2physeqhapmergedDate = merge_samples(HTACO2physeqhap, "Date")
SD= merge_samples(sample_data(HTACO2physeqhap), "Date")
print(SD[, "Date"])
print(HTACO2physeqhapmergedDate)
sample_names(HTACO2physeqhapmergedDate)
sample_data(HTACO2physeqhapmergedDate)$Date
identical(SD,sample_data(HTACO2physeqhapmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqhapmergedDate)$Date <- factor(sample_names(HTACO2physeqhapmergedDate))
#Getting abundance
hapHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqhapmergedDate), "matrix")
View(hapHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
hapHTACO2_abund_mergedDate_abs_pre<-hapHTACO2_abund_mergedDate
write.csv(hapHTACO2_abund_mergedDate, "hapHTACO2_abund_mergedDate.csv", row.names = F)

###Subset with Tethyida
HTACO2physeqteth=subset_taxa(HTACO2physeqpor,Order=="Tethyida")
HTACO2physeqtethmergedDate = merge_samples(HTACO2physeqteth, "Date")
SD= merge_samples(sample_data(HTACO2physeqteth), "Date")
print(SD[, "Date"])
print(HTACO2physeqtethmergedDate)
sample_names(HTACO2physeqtethmergedDate)
sample_data(HTACO2physeqtethmergedDate)$Date
identical(SD,sample_data(HTACO2physeqtethmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTACO2physeqtethmergedDate)$Date <- factor(sample_names(HTACO2physeqtethmergedDate))
#Getting abundance
tethHTACO2_abund_mergedDate=as(otu_table(HTACO2physeqtethmergedDate), "matrix")
View(tethHTACO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tethHTACO2_abund_mergedDate_abs_pre<-tethHTACO2_abund_mergedDate
write.csv(tethHTACO2_abund_mergedDate, "tethHTACO2_abund_mergedDate.csv", row.names = F)

#merge all of the csv per sponge group and bring in to R
poriferabyclasorder_HTACO2_abund_mergedDate<-read.csv("poriferabyclasorder_HTACO2_abund_mergedDate.csv") #import .csv poriferabyclasorder_Intake_abund_mergedDate.csv
View(poriferabyclasorder_HTACO2_abund_mergedDate)
poriferabyclasorder_HTACO2_abund_mergedDate[is.na(poriferabyclasorder_HTACO2_abund_mergedDate)] <- 0
View(poriferabyclasorder_HTACO2_abund_mergedDate)
Porif_HTACO2_class_order<-poriferabyclasorder_HTACO2_abund_mergedDate
names(Porif_HTACO2_class_order)
####rarefraction curve###
Porif_HTACO2_class_order_all<-Porif_HTACO2_class_order[2:94]
curve_Porif_HTACO2_class_order_all = specaccum(Porif_HTACO2_class_order_all, method = "rarefaction", 
                                               permutations = 100)
#subset each habitat into its own df

Porif_HTACO2_class_order%>% filter(Site == "Calcarea") -> Calcarea
names(Calcarea)
Porif_HTACO2_class_order%>% filter(Site == "Homoscleromorpha") -> Homoscleromorpha
Porif_HTACO2_class_order%>% filter(Site == "Keratosa") -> Keratosa
Porif_HTACO2_class_order%>% filter(Site == "Tetractinellida") -> Tetractinellida
Porif_HTACO2_class_order%>% filter(Site == "Suberitida") -> Suberitida
Porif_HTACO2_class_order%>% filter(Site == "Poecilosclerida") -> Poecilosclerida
Porif_HTACO2_class_order%>% filter(Site == "Haplosclerida") -> Haplosclerida
Porif_HTACO2_class_order%>% filter(Site == "Tethyida") -> Tethyida
Porif_HTACO2_class_order%>% filter(Site == "Dictyoceratida") -> Dictyoceratida
Porif_HTACO2_class_order%>% filter(Site == "Dendroceratida") -> Dendroceratida

#species accumulation curve for each habitat using all sponges
curve_HTACO2_calc = specaccum(Calcarea[, 2:94], method = "rarefaction")
curve_HTACO2_hom = specaccum(Homoscleromorpha[, 2:94], method = "rarefaction")
curve_HTACO2_ker = specaccum(Keratosa[, 2:94], method = "rarefaction")
curve_HTACO2_tetr = specaccum(Tetractinellida[, 2:94], method = "rarefaction")
curve_HTACO2_sub = specaccum(Suberitida[, 2:94], method = "rarefaction")
curve_HTACO2_poec = specaccum(Poecilosclerida[, 2:94], method = "rarefaction")
curve_HTACO2_hap = specaccum(Haplosclerida[, 2:94], method = "rarefaction")
curve_HTACO2_teth = specaccum(Tethyida[, 2:94], method = "rarefaction")
curve_HTACO2_Dicty = specaccum(Dictyoceratida[, 2:94], method = "rarefaction")
curve_HTACO2_Dendr = specaccum(Dendroceratida[, 2:94], method = "rarefaction")
#figure out color palet "display.brewer.pal(n = 8, name = 'RdBu')"
#got the error Error in plot.new() : figure margins too large" so I typed the following command "par(mar=c(5,5,3,1))"" # this fixed the margins and shows axis vallues. 
#plot curve_all first by dates
par(mar=c(5,5,3,1))
plot(curve_HTACO2_calc, xvar=c("individuals"),ylim=c(0,15),xlim=c(0,800),col="#B2182B", lwd=2, ci.lty=0, ci.col="#F4A582",
     main = "Default: Prettier CI")

#then plot the rest
plot(curve_HTACO2_hom, add = TRUE,xvar=c("individuals"), col="#2166AC", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI") #col is COLOUR setting, so change it to something else if you 

plot(curve_HTACO2_ker, add = TRUE,xvar=c("individuals"), col="#E1BE6A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTACO2_tetr, add = TRUE,xvar=c("individuals"), col="#40B0A6", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTACO2_sub, add = TRUE,xvar=c("individuals"), col="#E66100", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTACO2_poec, add = TRUE,xvar=c("individuals"), col="#5D3A9B", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTACO2_hap, add = TRUE,xvar=c("individuals"), col="#1AFF1A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTACO2_teth, add = TRUE,xvar=c("individuals"), col="#994F00", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")
legend(500, 1, legend=c("Calcareous", "Homoscleromorpha"),
       col=c("#B2182B", "#2166AC"), lty=1:1, cex=0.8)

#decided to exclude Dictyo Dendro and Verong since counts were minimal.
```

#rarefraction by TMT=ATHCO2

```{r setup, include=FALSE}
#Subset data by TMT="ATHCO2" 
ATHCO2physeqpor=subset_samples(physeqpor,TMT=="ATHCO2")
#OTHER MS For ARMS vs ReEf paper phylogenetic analysis etc. 
ATHCO2physeqpormergedDate = merge_samples(ATHCO2physeqpor, "Date")
SD= merge_samples(sample_data(ATHCO2physeqpor), "Date")
print(SD[, "Date"])
print(ATHCO2physeqpormergedDate)
sample_names(ATHCO2physeqpormergedDate)
sample_data(ATHCO2physeqpormergedDate)$Date
identical(SD,sample_data(ATHCO2physeqpormergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqpormergedDate)$Date <- factor(sample_names(ATHCO2physeqpormergedDate))
#Getting abundance
poriferaATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqpormergedDate), "matrix")
view(poriferaATHCO2_abund_mergedDate)
write.csv(poriferaATHCO2_abund_mergedDate, "poriferaATHCO2_abund_mergedDate.csv", row.names = F)

#incorporating abundance changing it to presence absence data
poriferaATHCO2_abund_mergedDate_abs_pre<-poriferaATHCO2_abund_mergedDate

poriferaATHCO2_abund_mergedDate_abs_pre[poriferaATHCO2_abund_mergedDate_abs_pre >0] <-1 #making a presence absence data frame
View(poriferaATHCO2_abund_mergedDate_abs_pre)
write.csv(poriferaATHCO2_abund_mergedDate_abs_pre, "poriferaATHCO2_abund_mergedDate_abs_pre.csv", row.names = F)

#Getting abundance
poriferaATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqpormergedDate), "matrix")
View(poriferaATHCO2_abund_mergedunit)
write.csv(poriferaATHCO2_abund_mergedunit, "poriferaATHCO2_abund_mergedunit.csv", row.names = F)

###subsetting samples by "ATHCO2"" and then merging sample by "Unit"" 
#which will combine the total diversity of the experiment for this treatment and calculate mean values of diversity#####
ATHCO2physeqpor=subset_samples(physeqpor,TMT=="ATHCO2")
#OTHER MS For ARMS vs ReEf paper phylogenetic analysis etc. 
ATHCO2physeqpormergedunit = merge_samples(ATHCO2physeqpor, "UNIT")
SD= merge_samples(sample_data(ATHCO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(ATHCO2physeqpormergedunit)
sample_names(ATHCO2physeqpormergedunit)
sample_data(ATHCO2physeqpormergedunit)$UNIT
identical(SD,sample_data(ATHCO2physeqpormergedunit))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqpormergedunit)$UNIT <- factor(sample_names(ATHCO2physeqpormergedunit))

#Determine diversity within Porifera ATHCO2.
Richness_ATHCO2physeqpormergedunit<-estimate_richness(ATHCO2physeqpormergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqpormergedunit, "Richness_ATHCO2physeqpormergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqpormergedunit)=sample_names(ATHCO2physeqpormergedunit)
Richness_ATHCO2physeqpormergedunit<-dfRowName(Richness_ATHCO2physeqpormergedunit, name ="UNIT")
View(Richness_ATHCO2physeqpormergedunit)
#need to add a column "TMT" to indicate "ATHCO2"
Richness_ATHCO2physeqpormergedunit<-cbind(Richness_ATHCO2physeqpormergedunit, TMT='ATHCO2')
View(Richness_ATHCO2physeqpormergedunit)
#Calculating mean values of Observed, SHannon, Simpson for the ATHCO2. which will be shown in table xx diversity. 
#Observed
Richness_ATHCO2physeqpormergedunit_sum<-summarySE(Richness_ATHCO2physeqpormergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqpormergedunit_sum)
#Shannon
Shannon_ATHCO2physeqpormergedunit_sum<-summarySE(Richness_ATHCO2physeqpormergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqpormergedunit_sum)
#Simpson
Simpson_ATHCO2physeqpormergedunit_sum<-summarySE(Richness_ATHCO2physeqpormergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqpormergedunit_sum)

#Getting abundance
poriferaATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqpormergedunit), "matrix")
View(poriferaATHCO2_abund_mergedunit)
write.csv(poriferaATHCO2_abund_mergedunit, "poriferaATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
poriferaATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(poriferaATHCO2_abund_mergedunit), data.frame(poriferaATHCO2_abund_mergedunit, row.names=NULL))
View(poriferaATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poriferaATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
poriferaATHCO2_abund_mergedunit_asdf<-as_data_frame(poriferaATHCO2_abund_mergedunit) 
View(poriferaATHCO2_abund_mergedunit_asdf)
#Sum rows
poriferaATHCO2_abund_mergedunit_asdf_total<-cbind(poriferaATHCO2_abund_mergedunit_asdf, total =rowSums(poriferaATHCO2_abund_mergedunit_asdf)) 
View(poriferaATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poriferaATHCO2_abund_mergedunit_rn_to_col" to "poriferaATHCO2_abund_mergedunit_rn_to_col_total"
poriferaATHCO2_abund_mergedunit_asdf_total$Unit<-poriferaATHCO2_abund_mergedunit_rn_to_col$Unit
names(poriferaATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 94. 
poriferaATHCO2_abund_mergedunit_asdf_total_clean<-poriferaATHCO2_abund_mergedunit_asdf_total[,-(1:93)]
View(poriferaATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
poriferaATHCO2_abund_mergedunit_asdf_total_clean<-cbind(poriferaATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
names(poriferaATHCO2_abund_mergedunit_asdf_total_clean)
View(poriferaATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(poriferaATHCO2_abund_mergedunit_asdf_total_clean, "poriferaATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poriferaATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poriferaATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poriferaATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Now calculate diversity within the different sponge subgroups within "ATHCO2physeqpormergedunit"" starting with Class Homoscleromorpha
ATHCO2physeqhomomergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Class=="Homoscleromorpha")
Richness_ATHCO2physeqhomomergedunit<-estimate_richness(ATHCO2physeqhomomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_ATHCO2physeqhomomergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqhomomergedunit, "Richness_ATHCO2physeqhomomergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqhomomergedunit)=sample_names(ATHCO2physeqhomomergedunit)
View(Richness_ATHCO2physeqhomomergedunit)
Richness_ATHCO2physeqhomomergedunit<-dfRowName(Richness_ATHCO2physeqhomomergedunit, name ="UNIT")
View(Richness_ATHCO2physeqhomomergedunit)
#Observed #add column with TMT name ATHCO2
Richness_ATHCO2physeqhomomergedunit <- Richness_ATHCO2physeqhomomergedunit %>%
  add_column(TMT = "ATHCO2")
Richness_ATHCO2physeqhomomergedunit
Richness_ATHCO2physeqhomomergedunit_sum<-summarySE(Richness_ATHCO2physeqhomomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqhomomergedunit_sum)
#Shannon
Shannon_ATHCO2physeqhomomergedunit_sum<-summarySE(Richness_ATHCO2physeqhomomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqhomomergedunit_sum)
#Simpson
Simpson_ATHCO2physeqhomomergedunit_sum<-summarySE(Richness_ATHCO2physeqhomomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqhomomergedunit_sum)

#Getting abundance
homoATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqhomomergedunit), "matrix")
View(homoATHCO2_abund_mergedunit)
write.csv(homoATHCO2_abund_mergedunit, "homoATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
homoATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(homoATHCO2_abund_mergedunit), data.frame(homoATHCO2_abund_mergedunit, row.names=NULL))
View(homoATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(homoATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
homoATHCO2_abund_mergedunit_asdf<-as_data_frame(homoATHCO2_abund_mergedunit) 
View(homoATHCO2_abund_mergedunit_asdf)
#Sum rows
homoATHCO2_abund_mergedunit_asdf_total<-cbind(homoATHCO2_abund_mergedunit_asdf, total =rowSums(homoATHCO2_abund_mergedunit_asdf)) 
View(homoATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "homoATHCO2_abund_mergedunit_rn_to_col" to "homoATHCO2_abund_mergedunit_rn_to_col_total"
homoATHCO2_abund_mergedunit_asdf_total$Unit<-homoATHCO2_abund_mergedunit_rn_to_col$Unit
View(homoATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
homoATHCO2_abund_mergedunit_asdf_total_clean<-homoATHCO2_abund_mergedunit_asdf_total[,-(1:8)]
View(homoATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
homoATHCO2_abund_mergedunit_asdf_total_clean<-cbind(homoATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(homoATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(homoATHCO2_abund_mergedunit_asdf_total_clean, "homoATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
homoATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(homoATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(homoATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Demospongiae
ATHCO2physeqdemmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Class=="Demospongiae")
Richness_ATHCO2physeqdemmergedunit<-estimate_richness(ATHCO2physeqdemmergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_ATHCO2physeqdemmergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqdemmergedunit, "Richness_ATHCO2physeqdemmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqdemmergedunit)=sample_names(ATHCO2physeqdemmergedunit)
View(Richness_ATHCO2physeqdemmergedunit)
Richness_ATHCO2physeqdemmergedunit<-dfRowName(Richness_ATHCO2physeqdemmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqdemmergedunit)
#Observed #add column with TMT name ATHCO2
Richness_ATHCO2physeqdemmergedunit <- Richness_ATHCO2physeqdemmergedunit %>%
  add_column(TMT = "ATHCO2")
Richness_ATHCO2physeqdemmergedunit
Richness_ATHCO2physeqdemmergedunit_sum<-summarySE(Richness_ATHCO2physeqdemmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqdemmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqdemmergedunit_sum<-summarySE(Richness_ATHCO2physeqdemmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqdemmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqdemmergedunit_sum<-summarySE(Richness_ATHCO2physeqdemmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqdemmergedunit_sum)

#Getting abundance
demATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqdemmergedunit), "matrix")
View(demATHCO2_abund_mergedunit)
write.csv(demATHCO2_abund_mergedunit, "demATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
demATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(demATHCO2_abund_mergedunit), data.frame(demATHCO2_abund_mergedunit, row.names=NULL))
View(demATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(demATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
demATHCO2_abund_mergedunit_asdf<-as_data_frame(demATHCO2_abund_mergedunit) 
View(demATHCO2_abund_mergedunit_asdf)
#Sum rows
demATHCO2_abund_mergedunit_asdf_total<-cbind(demATHCO2_abund_mergedunit_asdf, total =rowSums(demATHCO2_abund_mergedunit_asdf)) 
View(demATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "demATHCO2_abund_mergedunit_rn_to_col" to "demATHCO2_abund_mergedunit_rn_to_col_total"
demATHCO2_abund_mergedunit_asdf_total$Unit<-demATHCO2_abund_mergedunit_rn_to_col$Unit
View(demATHCO2_abund_mergedunit_asdf_total)
names(demATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 67. 
demATHCO2_abund_mergedunit_asdf_total_clean<-demATHCO2_abund_mergedunit_asdf_total[,-(1:66)]
names(demATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
demATHCO2_abund_mergedunit_asdf_total_clean<-cbind(demATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(demATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(demATHCO2_abund_mergedunit_asdf_total_clean, "demATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
demATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(demATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(demATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Dictyoceratida
#Now calculate diversity within the different sponge subgroups within "ATHCO2physeqpormergedunit"" starting with Order Dictyoceratida
ATHCO2physeqdictyomergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Dictyoceratida")
Richness_ATHCO2physeqdictyomergedunit<-estimate_richness(ATHCO2physeqdictyomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_ATHCO2physeqdictyomergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqdictyomergedunit, "Richness_ATHCO2physeqdictyomergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqdictyomergedunit)=sample_names(ATHCO2physeqdictyomergedunit)
View(Richness_ATHCO2physeqdictyomergedunit)
Richness_ATHCO2physeqdictyomergedunit<-dfRowName(Richness_ATHCO2physeqdictyomergedunit, name ="UNIT")
View(Richness_ATHCO2physeqdictyomergedunit)
#Observed #add column with TMT name ATHCO2
Richness_ATHCO2physeqdictyomergedunit <- Richness_ATHCO2physeqdictyomergedunit %>%
  add_column(TMT = "ATHCO2")
Richness_ATHCO2physeqdictyomergedunit
Richness_ATHCO2physeqdictyomergedunit_sum<-summarySE(Richness_ATHCO2physeqdictyomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqdictyomergedunit_sum)
#Shannon
Shannon_ATHCO2physeqdictyomergedunit_sum<-summarySE(Richness_ATHCO2physeqdictyomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqdictyomergedunit_sum)
#Simpson
Simpson_ATHCO2physeqdictyomergedunit_sum<-summarySE(Richness_ATHCO2physeqdictyomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqdictyomergedunit_sum)
#Getting abundance
dictyoATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqdictyomergedunit), "matrix")
View(dictyoATHCO2_abund_mergedunit)
write.csv(dictyoATHCO2_abund_mergedunit, "dictyoATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dictyoATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(dictyoATHCO2_abund_mergedunit), data.frame(dictyoATHCO2_abund_mergedunit, row.names=NULL))
View(dictyoATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dictyoATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dictyoATHCO2_abund_mergedunit_asdf<-as_data_frame(dictyoATHCO2_abund_mergedunit) 
View(dictyoATHCO2_abund_mergedunit_asdf)
#Sum rows
dictyoATHCO2_abund_mergedunit_asdf_total<-cbind(dictyoATHCO2_abund_mergedunit_asdf, total =rowSums(dictyoATHCO2_abund_mergedunit_asdf)) 
View(dictyoATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dictyoATHCO2_abund_mergedunit_rn_to_col" to "dictyoATHCO2_abund_mergedunit_rn_to_col_total"
dictyoATHCO2_abund_mergedunit_asdf_total$Unit<-dictyoATHCO2_abund_mergedunit_rn_to_col$Unit
View(dictyoATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dictyoATHCO2_abund_mergedunit_asdf_total_clean<-dictyoATHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dictyoATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
dictyoATHCO2_abund_mergedunit_asdf_total_clean<-cbind(dictyoATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(dictyoATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(dictyoATHCO2_abund_mergedunit_asdf_total_clean, "dictyoATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dictyoATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dictyoATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dictyoATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Dendroceratida
#Now calculate diversity within the different sponge subgroups within "ATHCO2physeqpormergedunit"" starting with Order Dendroceratida
ATHCO2physeqdendromergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Dendroceratida")
Richness_ATHCO2physeqdendromergedunit<-estimate_richness(ATHCO2physeqdendromergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_ATHCO2physeqdendromergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqdendromergedunit, "Richness_ATHCO2physeqdendromergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqdendromergedunit)=sample_names(ATHCO2physeqdendromergedunit)
View(Richness_ATHCO2physeqdendromergedunit)
Richness_ATHCO2physeqdendromergedunit<-dfRowName(Richness_ATHCO2physeqdendromergedunit, name ="UNIT")
View(Richness_ATHCO2physeqdendromergedunit)
#Observed #add column with TMT name ATHCO2
Richness_ATHCO2physeqdendromergedunit <- Richness_ATHCO2physeqdendromergedunit %>%
  add_column(TMT = "ATHCO2")
Richness_ATHCO2physeqdendromergedunit_sum<-summarySE(Richness_ATHCO2physeqdendromergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqdendromergedunit_sum)
#Shannon
Shannon_ATHCO2physeqdendromergedunit_sum<-summarySE(Richness_ATHCO2physeqdendromergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqdendromergedunit_sum)
#Simpson
Simpson_ATHCO2physeqdendromergedunit_sum<-summarySE(Richness_ATHCO2physeqdendromergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqdendromergedunit_sum)
#Getting abundance
dendroATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqdendromergedunit), "matrix")
View(dendroATHCO2_abund_mergedunit)
write.csv(dendroATHCO2_abund_mergedunit, "dendroATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dendroATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(dendroATHCO2_abund_mergedunit), data.frame(dendroATHCO2_abund_mergedunit, row.names=NULL))
View(dendroATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dendroATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dendroATHCO2_abund_mergedunit_asdf<-as_data_frame(dendroATHCO2_abund_mergedunit) 
View(dendroATHCO2_abund_mergedunit_asdf)
#Sum rows
dendroATHCO2_abund_mergedunit_asdf_total<-cbind(dendroATHCO2_abund_mergedunit_asdf, total =rowSums(dendroATHCO2_abund_mergedunit_asdf)) 
View(dendroATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dendroATHCO2_abund_mergedunit_rn_to_col" to "dendroATHCO2_abund_mergedunit_rn_to_col_total"
dendroATHCO2_abund_mergedunit_asdf_total$Unit<-dendroATHCO2_abund_mergedunit_rn_to_col$Unit
View(dendroATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dendroATHCO2_abund_mergedunit_asdf_total_clean<-dendroATHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dendroATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
dendroATHCO2_abund_mergedunit_asdf_total_clean<-cbind(dendroATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(dendroATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(dendroATHCO2_abund_mergedunit_asdf_total_clean, "dendroATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dendroATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dendroATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dendroATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Calcarea
ATHCO2physeqcalcmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Class=="Calcarea")
Richness_ATHCO2physeqcalcmergedunit<-estimate_richness(ATHCO2physeqcalcmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqcalcmergedunit, "Richness_ATHCO2physeqcalcmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqcalcmergedunit)=sample_names(ATHCO2physeqcalcmergedunit)
Richness_ATHCO2physeqcalcmergedunit<-dfRowName(Richness_ATHCO2physeqcalcmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqcalcmergedunit)
Richness_ATHCO2physeqcalcmergedunit<-cbind(Richness_ATHCO2physeqcalcmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqcalcmergedunit_sum<-summarySE(Richness_ATHCO2physeqcalcmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqcalcmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqcalcmergedunit_sum<-summarySE(Richness_ATHCO2physeqcalcmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqcalcmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqcalcmergedunit_sum<-summarySE(Richness_ATHCO2physeqcalcmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqcalcmergedunit_sum)
#Getting abundance
calcATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqcalcmergedunit), "matrix")
View(calcATHCO2_abund_mergedunit)
write.csv(calcATHCO2_abund_mergedunit, "calcATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
calcATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(calcATHCO2_abund_mergedunit), data.frame(calcATHCO2_abund_mergedunit, row.names=NULL))
View(calcATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(calcATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
calcATHCO2_abund_mergedunit_asdf<-as_data_frame(calcATHCO2_abund_mergedunit) 
View(calcATHCO2_abund_mergedunit_asdf)
#Sum rows
calcATHCO2_abund_mergedunit_asdf_total<-cbind(calcATHCO2_abund_mergedunit_asdf, total =rowSums(calcATHCO2_abund_mergedunit_asdf)) 
View(calcATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "calcATHCO2_abund_mergedunit_rn_to_col" to "calcATHCO2_abund_mergedunit_rn_to_col_total"
calcATHCO2_abund_mergedunit_asdf_total$Unit<-calcATHCO2_abund_mergedunit_rn_to_col$Unit
View(calcATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
calcATHCO2_abund_mergedunit_asdf_total_clean<-calcATHCO2_abund_mergedunit_asdf_total[,-(1:19)]
View(calcATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
calcATHCO2_abund_mergedunit_asdf_total_clean<-cbind(calcATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(calcATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(calcATHCO2_abund_mergedunit_asdf_total_clean, "calcATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
calcATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(calcATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(calcATHCO2_abund_mergedunit_asdf_total_clean_sum)


#Keratosa
ATHCO2physeqkermergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Subclass=="Keratosa")
Richness_ATHCO2physeqkermergedunit<-estimate_richness(ATHCO2physeqkermergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqkermergedunit, "Richness_ATHCO2physeqkermergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqkermergedunit)=sample_names(ATHCO2physeqkermergedunit)
Richness_ATHCO2physeqkermergedunit<-dfRowName(Richness_ATHCO2physeqkermergedunit, name ="UNIT")
View(Richness_ATHCO2physeqkermergedunit)
Richness_ATHCO2physeqkermergedunit<-cbind(Richness_ATHCO2physeqkermergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqkermergedunit_sum<-summarySE(Richness_ATHCO2physeqkermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqkermergedunit_sum)
#Shannon
Shannon_ATHCO2physeqkermergedunit_sum<-summarySE(Richness_ATHCO2physeqkermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqkermergedunit_sum)
#Simpson
Simpson_ATHCO2physeqkermergedunit_sum<-summarySE(Richness_ATHCO2physeqkermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqkermergedunit_sum)

#Getting abundance
kerATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqkermergedunit), "matrix")
View(kerATHCO2_abund_mergedunit)
write.csv(kerATHCO2_abund_mergedunit, "kerATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
kerATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(kerATHCO2_abund_mergedunit), data.frame(kerATHCO2_abund_mergedunit, row.names=NULL))
View(kerATHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(kerATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
kerATHCO2_abund_mergedunit_asdf<-as_data_frame(kerATHCO2_abund_mergedunit) 
View(kerATHCO2_abund_mergedunit_asdf)
#Sum rows
kerATHCO2_abund_mergedunit_asdf_total<-cbind(kerATHCO2_abund_mergedunit_asdf, total =rowSums(kerATHCO2_abund_mergedunit_asdf)) 
View(kerATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "kerATHCO2_abund_mergedunit_rn_to_col" to "kerATHCO2_abund_mergedunit_rn_to_col_total"
kerATHCO2_abund_mergedunit_asdf_total$Unit<-kerATHCO2_abund_mergedunit_rn_to_col$Unit
View(kerATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
kerATHCO2_abund_mergedunit_asdf_total_clean<-kerATHCO2_abund_mergedunit_asdf_total[,-(1:6)]
View(kerATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
kerATHCO2_abund_mergedunit_asdf_total_clean<-cbind(kerATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(kerATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(kerATHCO2_abund_mergedunit_asdf_total_clean, "kerATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
kerATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(kerATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(kerATHCO2_abund_mergedunit_asdf_total_clean_sum)

#abundance
keralltmtmerged_abun<-rbind(kerHTHCO2_abund_mergedunit_asdf_total_clean, kerATHCO2_abund_mergedunit_asdf_total_clean, kerAmbient_abund_mergedunit_asdf_total_clean,kerHTACO2_abund_mergedunit_asdf_total_clean)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 
keralltmtmerged_abun_df<-as.data.frame(keralltmtmerged_abun) #Had to change from matrix to dataframe
write_csv(keralltmtmerged_abun_df, "keralltmtmerged_abun_df.csv")
keralltmtmergeddfabun<-read.csv("keralltmtmerged_abun_df.csv") #import .csv
keralltmtmergeddfabun
#adding pco2 and temp as columns
keralltmtmergeddfabun <- keralltmtmergeddfabun %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
head(keralltmtmergeddfabun)
class(keralltmtmergeddfabun)
names(keralltmtmergeddfabun)
colnames(keralltmtmergeddfabun)[4] <- "temp" #change colnames of   TMT2 to temp
names(keralltmtmergeddfabun)
keralltmtmergeddfabun$temp <- revalue(keralltmtmergeddfabun$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(keralltmtmergeddfabun)
keralltmtmergeddfabun <- keralltmtmergeddfabun %>% 
  mutate(TMT2 = TMT)
names(keralltmtmergeddfabun)
colnames(keralltmtmergeddfabun)[5] <- "pco2" #change colnames of   TMT2 to temp
names(keralltmtmergeddfabun)
keralltmtmergeddfabun$pco2 <- revalue(keralltmtmergeddfabun$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))

#violin plot
#violin plots guided by https://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
#with dots
#effect of temp
keralltempmergeddabuntemp_viol3<-ggplot(keralltmtmergeddfabun, aes(x = temp, y = total, fill = temp)) + scale_x_discrete(limits=c("Present", "Future"))+theme_classic()+geom_violin(alpha = 0.5) +  ylim(0,320)+
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(y = "Abundance", x = element_blank())+scale_fill_manual("legend", values=c("Future" = "grey67","Present" ="cornsilk2"))

keralltempmergeddabuntemp_viol3
#with quartile and mean 
keralltempmergeddfabuntemp_viol4<-keralltempmergeddabuntemp_viol3 + geom_boxplot(width=0.1)
keralltempmergeddfabuntemp_viol4
keralltempmergeddabuntemp_viol3_tempsum<-summarySE(keralltmtmergeddfabun, measurevar = "total", groupvars=c("temp"))
keralltempmergeddabuntemp_viol3_tempsum


#effect of pCO2
kerallpco2mergeddabunpco2_viol3<-ggplot(keralltmtmergeddfabun, aes(x = pco2, y = total, fill = pco2)) + scale_x_discrete(limits=c("Present", "Future"))+theme_classic()+geom_violin(alpha = 0.5)+  ylim(0,320)+
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(y = "Abundance", x = element_blank())+scale_fill_manual("legend", values=c("Future" = "grey67","Present" ="cornsilk2"))
kerallpco2mergeddabunpco2_viol3
#with quartile and mean 
kerallpco2mergeddfabunpco2_viol4<-kerallpco2mergeddabunpco2_viol3 + geom_boxplot(width=0.1)
kerallpco2mergeddfabunpco2_viol4
keralltempmergeddabuntemp_viol3_pco2sum<-summarySE(keralltmtmergeddfabun, measurevar = "total", groupvars=c("pco2"))
keralltempmergeddabuntemp_viol3_pco2sum
#Subclass==Verongimorpha
ATHCO2physeqvermergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Subclass=="Verongimorpha")
Richness_ATHCO2physeqvermergedunit<-estimate_richness(ATHCO2physeqvermergedunit, measures = c("Observed", "Shannon", "Simpson"))
write.csv(Richness_ATHCO2physeqvermergedunit, "Richness_ATHCO2physeqvermergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqvermergedunit)=sample_names(ATHCO2physeqvermergedunit)
Richness_ATHCO2physeqvermergedunit<-dfRowName(Richness_ATHCO2physeqvermergedunit, name ="UNIT")
View(Richness_ATHCO2physeqvermergedunit)
Richness_ATHCO2physeqvermergedunit<-cbind(Richness_ATHCO2physeqvermergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqvermergedunit_sum<-summarySE(Richness_ATHCO2physeqvermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqvermergedunit_sum)
#Shannon
Shannon_ATHCO2physeqvermergedunit_sum<-summarySE(Richness_ATHCO2physeqvermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqvermergedunit_sum)
#Simpson
Simpson_ATHCO2physeqvermergedunit_sum<-summarySE(Richness_ATHCO2physeqvermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqvermergedunit_sum)
#Getting abundance
verATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqvermergedunit), "matrix")
View(verATHCO2_abund_mergedunit)
write.csv(verATHCO2_abund_mergedunit, "verATHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
verATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(verATHCO2_abund_mergedunit), data.frame(verATHCO2_abund_mergedunit, row.names=NULL))
View(verATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(verATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
verATHCO2_abund_mergedunit_asdf<-as_data_frame(verATHCO2_abund_mergedunit) 
View(verATHCO2_abund_mergedunit_asdf)
#Sum rows
verATHCO2_abund_mergedunit_asdf_total<-cbind(verATHCO2_abund_mergedunit_asdf, total =rowSums(verATHCO2_abund_mergedunit_asdf)) 
View(verATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "verATHCO2_abund_mergedunit_rn_to_col" to "verATHCO2_abund_mergedunit_rn_to_col_total"
verATHCO2_abund_mergedunit_asdf_total$Unit<-verATHCO2_abund_mergedunit_rn_to_col$Unit
View(verATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
verATHCO2_abund_mergedunit_asdf_total_clean<-verATHCO2_abund_mergedunit_asdf_total[,-(1:2)]
View(verATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
verATHCO2_abund_mergedunit_asdf_total_clean<-cbind(verATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(verATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(verATHCO2_abund_mergedunit_asdf_total_clean, "verATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
verATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(verATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(verATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Suberitida
ATHCO2physeqsubmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Suberitida")
Richness_ATHCO2physeqsubmergedunit<-estimate_richness(ATHCO2physeqsubmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqsubmergedunit, "Richness_ATHCO2physeqsubmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqsubmergedunit)=sample_names(ATHCO2physeqsubmergedunit)
Richness_ATHCO2physeqsubmergedunit<-dfRowName(Richness_ATHCO2physeqsubmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqsubmergedunit)
Richness_ATHCO2physeqsubmergedunit<-cbind(Richness_ATHCO2physeqsubmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqsubmergedunit_sum<-summarySE(Richness_ATHCO2physeqsubmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqsubmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqsubmergedunit_sum<-summarySE(Richness_ATHCO2physeqsubmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqsubmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqsubmergedunit_sum<-summarySE(Richness_ATHCO2physeqsubmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqsubmergedunit_sum)
#Getting abundance
subATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqsubmergedunit), "matrix")
View(subATHCO2_abund_mergedunit)
write.csv(subATHCO2_abund_mergedunit, "subATHCO2_abund_mergedunit.csv", row.names = F)
#Consubt Rowname into first column 
subATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(subATHCO2_abund_mergedunit), data.frame(subATHCO2_abund_mergedunit, row.names=NULL))
View(subATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(subATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#consubt to dataframe
subATHCO2_abund_mergedunit_asdf<-as_data_frame(subATHCO2_abund_mergedunit) 
View(subATHCO2_abund_mergedunit_asdf)
#Sum rows
subATHCO2_abund_mergedunit_asdf_total<-cbind(subATHCO2_abund_mergedunit_asdf, total =rowSums(subATHCO2_abund_mergedunit_asdf)) 
View(subATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "subATHCO2_abund_mergedunit_rn_to_col" to "subATHCO2_abund_mergedunit_rn_to_col_total"
subATHCO2_abund_mergedunit_asdf_total$Unit<-subATHCO2_abund_mergedunit_rn_to_col$Unit
View(subATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
subATHCO2_abund_mergedunit_asdf_total_clean<-subATHCO2_abund_mergedunit_asdf_total[,-(1:15)]
View(subATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
subATHCO2_abund_mergedunit_asdf_total_clean<-cbind(subATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(subATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(subATHCO2_abund_mergedunit_asdf_total_clean, "subATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
subATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(subATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(subATHCO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Tetractinellida
ATHCO2physeqtetrmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Tetractinellida")
Richness_ATHCO2physeqtetrmergedunit<-estimate_richness(ATHCO2physeqtetrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqtetrmergedunit, "Richness_ATHCO2physeqtetrmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqtetrmergedunit)=sample_names(ATHCO2physeqtetrmergedunit)
Richness_ATHCO2physeqtetrmergedunit<-dfRowName(Richness_ATHCO2physeqtetrmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqtetrmergedunit)
Richness_ATHCO2physeqtetrmergedunit<-cbind(Richness_ATHCO2physeqtetrmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqtetrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtetrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqtetrmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqtetrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtetrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqtetrmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqtetrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtetrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqtetrmergedunit_sum)
#Getting abundance
tetrATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqtetrmergedunit), "matrix")
View(tetrATHCO2_abund_mergedunit)
write.csv(tetrATHCO2_abund_mergedunit, "tetrATHCO2_abund_mergedunit.csv", row.names = F)
#Contetrt Rowname into first column 
tetrATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(tetrATHCO2_abund_mergedunit), data.frame(tetrATHCO2_abund_mergedunit, row.names=NULL))
View(tetrATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tetrATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contetrt to dataframe
tetrATHCO2_abund_mergedunit_asdf<-as_data_frame(tetrATHCO2_abund_mergedunit) 
View(tetrATHCO2_abund_mergedunit_asdf)
#Sum rows
tetrATHCO2_abund_mergedunit_asdf_total<-cbind(tetrATHCO2_abund_mergedunit_asdf, total =rowSums(tetrATHCO2_abund_mergedunit_asdf)) 
View(tetrATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tetrATHCO2_abund_mergedunit_rn_to_col" to "tetrATHCO2_abund_mergedunit_rn_to_col_total"
tetrATHCO2_abund_mergedunit_asdf_total$Unit<-tetrATHCO2_abund_mergedunit_rn_to_col$Unit
View(tetrATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
tetrATHCO2_abund_mergedunit_asdf_total_clean<-tetrATHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(tetrATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
tetrATHCO2_abund_mergedunit_asdf_total_clean<-cbind(tetrATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(tetrATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(tetrATHCO2_abund_mergedunit_asdf_total_clean, "tetrATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tetrATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tetrATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tetrATHCO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Poecilosclerida
ATHCO2physeqpoecrmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Poecilosclerida")
Richness_ATHCO2physeqpoecrmergedunit<-estimate_richness(ATHCO2physeqpoecrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqpoecrmergedunit, "Richness_ATHCO2physeqpoecrmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqpoecrmergedunit)=sample_names(ATHCO2physeqpoecrmergedunit)
Richness_ATHCO2physeqpoecrmergedunit<-dfRowName(Richness_ATHCO2physeqpoecrmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqpoecrmergedunit)
Richness_ATHCO2physeqpoecrmergedunit<-cbind(Richness_ATHCO2physeqpoecrmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_ATHCO2physeqpoecrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqpoecrmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_ATHCO2physeqpoecrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqpoecrmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_ATHCO2physeqpoecrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqpoecrmergedunit_sum)
#Getting abundance
poecrATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqpoecrmergedunit), "matrix")
View(poecrATHCO2_abund_mergedunit)
write.csv(poecrATHCO2_abund_mergedunit, "poecrATHCO2_abund_mergedunit.csv", row.names = F)
#Conpoecrt Rowname into first column 
poecrATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(poecrATHCO2_abund_mergedunit), data.frame(poecrATHCO2_abund_mergedunit, row.names=NULL))
View(poecrATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poecrATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conpoecrt to dataframe
poecrATHCO2_abund_mergedunit_asdf<-as_data_frame(poecrATHCO2_abund_mergedunit) 
View(poecrATHCO2_abund_mergedunit_asdf)
#Sum rows
poecrATHCO2_abund_mergedunit_asdf_total<-cbind(poecrATHCO2_abund_mergedunit_asdf, total =rowSums(poecrATHCO2_abund_mergedunit_asdf)) 
View(poecrATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poecrATHCO2_abund_mergedunit_rn_to_col" to "poecrATHCO2_abund_mergedunit_rn_to_col_total"
poecrATHCO2_abund_mergedunit_asdf_total$Unit<-poecrATHCO2_abund_mergedunit_rn_to_col$Unit
View(poecrATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
poecrATHCO2_abund_mergedunit_asdf_total_clean<-poecrATHCO2_abund_mergedunit_asdf_total[,-(1:12)]
View(poecrATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
poecrATHCO2_abund_mergedunit_asdf_total_clean<-cbind(poecrATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(poecrATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(poecrATHCO2_abund_mergedunit_asdf_total_clean, "poecrATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poecrATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poecrATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poecrATHCO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Tethyida
ATHCO2physeqtethrmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Tethyida")
Richness_ATHCO2physeqtethrmergedunit<-estimate_richness(ATHCO2physeqtethrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqtethrmergedunit, "Richness_ATHCO2physeqtethrmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqtethrmergedunit)=sample_names(ATHCO2physeqtethrmergedunit)
Richness_ATHCO2physeqtethrmergedunit<-dfRowName(Richness_ATHCO2physeqtethrmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqtethrmergedunit)
Richness_ATHCO2physeqtethrmergedunit<-cbind(Richness_ATHCO2physeqtethrmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqtethrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtethrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqtethrmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqtethrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtethrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqtethrmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqtethrmergedunit_sum<-summarySE(Richness_ATHCO2physeqtethrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqtethrmergedunit_sum)
#Getting abundance
tethrATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqtethrmergedunit), "matrix")
View(tethrATHCO2_abund_mergedunit)
write.csv(tethrATHCO2_abund_mergedunit, "tethrATHCO2_abund_mergedunit.csv", row.names = F)
#Contethrt Rowname into first column 
tethrATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(tethrATHCO2_abund_mergedunit), data.frame(tethrATHCO2_abund_mergedunit, row.names=NULL))
View(tethrATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tethrATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contethrt to dataframe
tethrATHCO2_abund_mergedunit_asdf<-as_data_frame(tethrATHCO2_abund_mergedunit) 
View(tethrATHCO2_abund_mergedunit_asdf)
#Sum rows
tethrATHCO2_abund_mergedunit_asdf_total<-cbind(tethrATHCO2_abund_mergedunit_asdf, total =rowSums(tethrATHCO2_abund_mergedunit_asdf)) 
View(tethrATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tethrATHCO2_abund_mergedunit_rn_to_col" to "tethrATHCO2_abund_mergedunit_rn_to_col_total"
tethrATHCO2_abund_mergedunit_asdf_total$Unit<-tethrATHCO2_abund_mergedunit_rn_to_col$Unit
View(tethrATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
tethrATHCO2_abund_mergedunit_asdf_total_clean<-tethrATHCO2_abund_mergedunit_asdf_total[,-(1:6)]
View(tethrATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
tethrATHCO2_abund_mergedunit_asdf_total_clean<-cbind(tethrATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(tethrATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(tethrATHCO2_abund_mergedunit_asdf_total_clean, "tethrATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tethrATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tethrATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tethrATHCO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Haplosclerida
ATHCO2physeqhaprmergedunit<-subset_taxa(ATHCO2physeqpormergedunit, Order=="Haplosclerida")
Richness_ATHCO2physeqhaprmergedunit<-estimate_richness(ATHCO2physeqhaprmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_ATHCO2physeqhaprmergedunit, "Richness_ATHCO2physeqhaprmergedunit.csv", row.names = F)
row.names(Richness_ATHCO2physeqhaprmergedunit)=sample_names(ATHCO2physeqhaprmergedunit)
Richness_ATHCO2physeqhaprmergedunit<-dfRowName(Richness_ATHCO2physeqhaprmergedunit, name ="UNIT")
View(Richness_ATHCO2physeqhaprmergedunit)
Richness_ATHCO2physeqhaprmergedunit<-cbind(Richness_ATHCO2physeqhaprmergedunit, TMT='ATHCO2')
#Observed
Richness_ATHCO2physeqhaprmergedunit_sum<-summarySE(Richness_ATHCO2physeqhaprmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_ATHCO2physeqhaprmergedunit_sum)
#Shannon
Shannon_ATHCO2physeqhaprmergedunit_sum<-summarySE(Richness_ATHCO2physeqhaprmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_ATHCO2physeqhaprmergedunit_sum)
#Simpson
Simpson_ATHCO2physeqhaprmergedunit_sum<-summarySE(Richness_ATHCO2physeqhaprmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_ATHCO2physeqhaprmergedunit_sum)
#Getting abundance
haprATHCO2_abund_mergedunit=as(otu_table(ATHCO2physeqhaprmergedunit), "matrix")
View(haprATHCO2_abund_mergedunit)
write.csv(haprATHCO2_abund_mergedunit, "haprATHCO2_abund_mergedunit.csv", row.names = F)
#Conhaprt Rowname into first column 
haprATHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(haprATHCO2_abund_mergedunit), data.frame(haprATHCO2_abund_mergedunit, row.names=NULL))
View(haprATHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(haprATHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conhaprt to dataframe
haprATHCO2_abund_mergedunit_asdf<-as_data_frame(haprATHCO2_abund_mergedunit) 
View(haprATHCO2_abund_mergedunit_asdf)
#Sum rows
haprATHCO2_abund_mergedunit_asdf_total<-cbind(haprATHCO2_abund_mergedunit_asdf, total =rowSums(haprATHCO2_abund_mergedunit_asdf)) 
View(haprATHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "haprATHCO2_abund_mergedunit_rn_to_col" to "haprATHCO2_abund_mergedunit_rn_to_col_total"
haprATHCO2_abund_mergedunit_asdf_total$Unit<-haprATHCO2_abund_mergedunit_rn_to_col$Unit
names(haprATHCO2_abund_mergedunit_asdf_total)
View(haprATHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 23. 
haprATHCO2_abund_mergedunit_asdf_total_clean<-haprATHCO2_abund_mergedunit_asdf_total[,-(1:22)]
View(haprATHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "ATHCO2"
haprATHCO2_abund_mergedunit_asdf_total_clean<-cbind(haprATHCO2_abund_mergedunit_asdf_total_clean, TMT='ATHCO2')
View(haprATHCO2_abund_mergedunit_asdf_total_clean)
write.csv(haprATHCO2_abund_mergedunit_asdf_total_clean, "haprATHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
haprATHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(haprATHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(haprATHCO2_abund_mergedunit_asdf_total_clean_sum)


#Rarefraction curve preliminary steps:Subset data within "ATHCO2DATE"to obtain the different subsponge groups.
#remember that when subsetting taxonomic assignments you must use "subset_taxa" starting with 
#Homoscleromorpha
ATHCO2physeqhomo=subset_taxa(ATHCO2physeqpor,Class=="Homoscleromorpha")
ATHCO2physeqhomomergedDate = merge_samples(ATHCO2physeqhomo, "Date")
SD= merge_samples(sample_data(ATHCO2physeqhomo), "Date")
print(SD[, "Date"])
print(ATHCO2physeqhomomergedDate)
sample_names(ATHCO2physeqhomomergedDate)
sample_data(ATHCO2physeqhomomergedDate)$Date
identical(SD,sample_data(ATHCO2physeqhomomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqhomomergedDate)$Date <- factor(sample_names(ATHCO2physeqhomomergedDate))
#Getting abundance
homoATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqhomomergedDate), "matrix")
View(homoATHCO2_abund_mergedDate)
#adding column at the end that shows TMT=ATHCO2
homoATHCO2_abund_mergedDate<-cbind(homoATHCO2_abund_mergedDate, TMT='ATHCO2')
head(homoATHCO2_abund_mergedDate)
write.csv(homoATHCO2_abund_mergedDate, "homoATHCO2_abund_mergedDate.csv", row.names = F)

#Demospongiae
ATHCO2physeqdem=subset_taxa(ATHCO2physeqpor,Class=="Demospongiae")
ATHCO2physeqdemmergedDate = merge_samples(ATHCO2physeqdem, "Date")
SD= merge_samples(sample_data(ATHCO2physeqdem), "Date")
print(SD[, "Date"])
print(ATHCO2physeqdemmergedDate)
sample_names(ATHCO2physeqdemmergedDate)
sample_data(ATHCO2physeqdemmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqdemmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqdemmergedDate)$Date <- factor(sample_names(ATHCO2physeqdemmergedDate))
#Getting abundance
demATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqdemmergedDate), "matrix")
View(demATHCO2_abund_mergedDate)
#adding column at the end that shows TMT=ATHCO2
demATHCO2_abund_mergedDate<-cbind(demATHCO2_abund_mergedDate, TMT='ATHCO2')
write.csv(demATHCO2_abund_mergedDate, "demATHCO2_abund_mergedDate.csv", row.names = F)

#Dictyoceratida
ATHCO2physeqdictyo=subset_taxa(ATHCO2physeqpor,Order=="Dictyoceratida")
ATHCO2physeqdictyomergedDate = merge_samples(ATHCO2physeqdictyo, "Date")
SD= merge_samples(sample_data(ATHCO2physeqdictyo), "Date")
print(SD[, "Date"])
print(ATHCO2physeqdictyomergedDate)
sample_names(ATHCO2physeqdictyomergedDate)
sample_data(ATHCO2physeqdictyomergedDate)$Date
identical(SD,sample_data(ATHCO2physeqdictyomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqdictyomergedDate)$Date <- factor(sample_names(ATHCO2physeqdictyomergedDate))
#Getting abundance
dictyoATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqdictyomergedDate), "matrix")
View(dictyoATHCO2_abund_mergedDate)
write.csv(dictyoATHCO2_abund_mergedDate, "dictyoATHCO2_abund_mergedDate.csv", row.names = F)

#Dendroceratida
ATHCO2physeqdendro=subset_taxa(ATHCO2physeqpor,Order=="Dendroceratida")
ATHCO2physeqdendromergedDate = merge_samples(ATHCO2physeqdendro, "Date")
SD= merge_samples(sample_data(ATHCO2physeqdendro), "Date")
print(SD[, "Date"])
print(ATHCO2physeqdendromergedDate)
sample_names(ATHCO2physeqdendromergedDate)
sample_data(ATHCO2physeqdendromergedDate)$Date
identical(SD,sample_data(ATHCO2physeqdendromergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqdendromergedDate)$Date <- factor(sample_names(ATHCO2physeqdendromergedDate))
#Getting abundance
dendroATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqdendromergedDate), "matrix")
View(dendroATHCO2_abund_mergedDate)
write.csv(dendroATHCO2_abund_mergedDate, "dendroATHCO2_abund_mergedDate.csv", row.names = F)


###Subset with Clathrinida 
ATHCO2physeqclat=subset_taxa(ATHCO2physeqpor,Order=="Clathrinida")
ATHCO2physeqclatrmergedDate = merge_samples(ATHCO2physeqclat, "Date")
SD= merge_samples(sample_data(ATHCO2physeqclat), "Date")
print(SD[, "Date"])
print(ATHCO2physeqclatrmergedDate)
sample_names(ATHCO2physeqclatrmergedDate)
sample_data(ATHCO2physeqclatrmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqclatrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqclatrmergedDate)$Date <- factor(sample_names(ATHCO2physeqclatrmergedDate))
#Getting abundance
clatATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqclatrmergedDate), "matrix")
View(clatATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
clatATHCO2_abund_mergedDate_abs_pre<-clatATHCO2_abund_mergedDate
write.csv(clatATHCO2_abund_mergedDate, "clatATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Leucosolenida
ATHCO2physeqleuc=subset_taxa(ATHCO2physeqpor,Order=="Leucosolenida")
ATHCO2physeqleucmergedDate = merge_samples(ATHCO2physeqleuc, "Date")
SD= merge_samples(sample_data(ATHCO2physeqleuc), "Date")
print(SD[, "Date"])
print(ATHCO2physeqleucmergedDate)
sample_names(ATHCO2physeqleucmergedDate)
sample_data(ATHCO2physeqleucmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqleucmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqleucmergedDate)$Date <- factor(sample_names(ATHCO2physeqleucmergedDate))
#Getting abundance
leucATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqleucmergedDate), "matrix")
View(leucATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
leucATHCO2_abund_mergedDate_abs_pre<-leucATHCO2_abund_mergedDate
write.csv(leucATHCO2_abund_mergedDate, "leucATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Class: Calcarea
ATHCO2physeqcalc=subset_taxa(ATHCO2physeqpor,Class=="Calcarea")
ATHCO2physeqcalcmergedDate = merge_samples(ATHCO2physeqcalc, "Date")
SD= merge_samples(sample_data(ATHCO2physeqcalc), "Date")
print(SD[, "Date"])
print(ATHCO2physeqcalcmergedDate)
sample_names(ATHCO2physeqcalcmergedDate)
sample_data(ATHCO2physeqcalcmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqcalcmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqcalcmergedDate)$Date <- factor(sample_names(ATHCO2physeqcalcmergedDate))
#Getting abundance
calcATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqcalcmergedDate), "matrix")
View(calcATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
calcATHCO2_abund_mergedDate_abs_pre<-calcATHCO2_abund_mergedDate
#adding column at the end that shows class=Demospongiae
calcATHCO2_abund_mergedDate<-cbind(calcATHCO2_abund_mergedDate, TMT='ATHCO2')
write.csv(calcATHCO2_abund_mergedDate, "calcATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Keratosa 
ATHCO2physeqker=subset_taxa(ATHCO2physeqpor,Subclass=="Keratosa")
ATHCO2physeqkermergedDate = merge_samples(ATHCO2physeqker, "Date")
SD= merge_samples(sample_data(ATHCO2physeqker), "Date")
print(SD[, "Date"])
print(ATHCO2physeqkermergedDate)
sample_names(ATHCO2physeqkermergedDate)
sample_data(ATHCO2physeqkermergedDate)$Date
identical(SD,sample_data(ATHCO2physeqkermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqkermergedDate)$Date <- factor(sample_names(ATHCO2physeqkermergedDate))
#Getting abundance
kerATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqkermergedDate), "matrix")
View(kerATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
kerATHCO2_abund_mergedDate_abs_pre<-kerATHCO2_abund_mergedDate
write.csv(kerATHCO2_abund_mergedDate, "kerATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Tetractinellida 
ATHCO2physeqtetr=subset_taxa(ATHCO2physeqpor,Order=="Tetractinellida")
ATHCO2physeqtetrmergedDate = merge_samples(ATHCO2physeqtetr, "Date")
SD= merge_samples(sample_data(ATHCO2physeqtetr), "Date")
print(SD[, "Date"])
print(ATHCO2physeqtetrmergedDate)
sample_names(ATHCO2physeqtetrmergedDate)
sample_data(ATHCO2physeqtetrmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqtetrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqtetrmergedDate)$Date <- factor(sample_names(ATHCO2physeqtetrmergedDate))
#Getting abundance
tetrATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqtetrmergedDate), "matrix")
View(tetrATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tetrATHCO2_abund_mergedDate_abs_pre<-tetrATHCO2_abund_mergedDate
write.csv(tetrATHCO2_abund_mergedDate, "tetrATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Suberitida 
ATHCO2physeqsub=subset_taxa(ATHCO2physeqpor,Order=="Suberitida")
ATHCO2physeqsubmergedDate = merge_samples(ATHCO2physeqsub, "Date")
SD= merge_samples(sample_data(ATHCO2physeqsub), "Date")
print(SD[, "Date"])
print(ATHCO2physeqsubmergedDate)
sample_names(ATHCO2physeqsubmergedDate)
sample_data(ATHCO2physeqsubmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqsubmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqsubmergedDate)$Date <- factor(sample_names(ATHCO2physeqsubmergedDate))
#Getting abundance
subATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqsubmergedDate), "matrix")
View(subATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
subATHCO2_abund_mergedDate_abs_pre<-subATHCO2_abund_mergedDate
write.csv(subATHCO2_abund_mergedDate, "subATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Poecilosclerida
ATHCO2physeqpoec=subset_taxa(ATHCO2physeqpor,Order=="Poecilosclerida")
ATHCO2physeqpoecmergedDate = merge_samples(ATHCO2physeqpoec, "Date")
SD= merge_samples(sample_data(ATHCO2physeqpoec), "Date")
print(SD[, "Date"])
print(ATHCO2physeqpoecmergedDate)
sample_names(ATHCO2physeqpoecmergedDate)
sample_data(ATHCO2physeqpoecmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqpoecmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqpoecmergedDate)$Date <- factor(sample_names(ATHCO2physeqpoecmergedDate))
#Getting abundance
poecATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqpoecmergedDate), "matrix")
View(poecATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poecATHCO2_abund_mergedDate_abs_pre<-poecATHCO2_abund_mergedDate
write.csv(poecATHCO2_abund_mergedDate, "poecATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Haplosclerida
ATHCO2physeqhap=subset_taxa(ATHCO2physeqpor,Order=="Haplosclerida")
ATHCO2physeqhapmergedDate = merge_samples(ATHCO2physeqhap, "Date")
SD= merge_samples(sample_data(ATHCO2physeqhap), "Date")
print(SD[, "Date"])
print(ATHCO2physeqhapmergedDate)
sample_names(ATHCO2physeqhapmergedDate)
sample_data(ATHCO2physeqhapmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqhapmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqhapmergedDate)$Date <- factor(sample_names(ATHCO2physeqhapmergedDate))
#Getting abundance
hapATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqhapmergedDate), "matrix")
View(hapATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
hapATHCO2_abund_mergedDate_abs_pre<-hapATHCO2_abund_mergedDate
write.csv(hapATHCO2_abund_mergedDate, "hapATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Tethyida
ATHCO2physeqteth=subset_taxa(ATHCO2physeqpor,Order=="Tethyida")
ATHCO2physeqtethmergedDate = merge_samples(ATHCO2physeqteth, "Date")
SD= merge_samples(sample_data(ATHCO2physeqteth), "Date")
print(SD[, "Date"])
print(ATHCO2physeqtethmergedDate)
sample_names(ATHCO2physeqtethmergedDate)
sample_data(ATHCO2physeqtethmergedDate)$Date
identical(SD,sample_data(ATHCO2physeqtethmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqtethmergedDate)$Date <- factor(sample_names(ATHCO2physeqtethmergedDate))
#Getting abundance
tethATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqtethmergedDate), "matrix")
View(tethATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tethATHCO2_abund_mergedDate_abs_pre<-tethATHCO2_abund_mergedDate
write.csv(tethATHCO2_abund_mergedDate, "tethATHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Verongida
ATHCO2physeqver=subset_taxa(ATHCO2physeqpor,Subclass=="Verongimorpha")
ATHCO2physeqvermergedDate = merge_samples(ATHCO2physeqver, "Date")
SD= merge_samples(sample_data(ATHCO2physeqver), "Date")
print(SD[, "Date"])
print(ATHCO2physeqvermergedDate)
sample_names(ATHCO2physeqvermergedDate)
sample_data(ATHCO2physeqvermergedDate)$Date
identical(SD,sample_data(ATHCO2physeqvermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(ATHCO2physeqvermergedDate)$Date <- factor(sample_names(ATHCO2physeqvermergedDate))
#Getting abundance
verATHCO2_abund_mergedDate=as(otu_table(ATHCO2physeqvermergedDate), "matrix")
View(verATHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
verATHCO2_abund_mergedDate_abs_pre<-verATHCO2_abund_mergedDate
write.csv(verATHCO2_abund_mergedDate, "verATHCO2_abund_mergedDate.csv", row.names = F)


#merge all of the csv per sponge group and bring in to R
poriferabyclasorder_ATHCO2_abund_mergedDate<-read.csv("poriferabyclasorder_ATHCO2_abund_mergedDate.csv") #import .csv poriferabyclasorder_Intake_abund_mergedDate.csv
names(poriferabyclasorder_ATHCO2_abund_mergedDate)
View(poriferabyclasorder_ATHCO2_abund_mergedDate)
poriferabyclasorder_ATHCO2_abund_mergedDate[is.na(poriferabyclasorder_ATHCO2_abund_mergedDate)] <- 0
View(poriferabyclasorder_ATHCO2_abund_mergedDate)
Porif_ATHCO2_class_order<-poriferabyclasorder_ATHCO2_abund_mergedDate
View(Porif_ATHCO2_class_order)
####rarefraction curve###
Porif_ATHCO2_class_order_all<-Porif_ATHCO2_class_order[2:94]
curve_Porif_ATHCO2_class_order_all = specaccum(Porif_ATHCO2_class_order_all, method = "rarefaction", 
                                               permutations = 100)
#subset each habitat into its own df
Porif_ATHCO2_class_order%>% filter(Site == "Calcarea") -> Calcarea
names(Calcarea)
Porif_ATHCO2_class_order%>% filter(Site == "Homoscleromorpha") -> Homoscleromorpha
Porif_ATHCO2_class_order%>% filter(Site == "Keratosa") -> Keratosa
Porif_ATHCO2_class_order%>% filter(Site == "Tetractinellida") -> Tetractinellida
Porif_ATHCO2_class_order%>% filter(Site == "Suberitida") -> Suberitida
Porif_ATHCO2_class_order%>% filter(Site == "Poecilosclerida") -> Poecilosclerida
Porif_ATHCO2_class_order%>% filter(Site == "Haplosclerida") -> Haplosclerida
Porif_ATHCO2_class_order%>% filter(Site == "Tethyida") -> Tethyida
Porif_ATHCO2_class_order%>% filter(Site == "Dictyoceratida") -> Dictyoceratida
Porif_ATHCO2_class_order%>% filter(Site == "Dendroceratida") -> Dendroceratida
Porif_ATHCO2_class_order%>% filter(Site == "Verongida") -> Verongida
#species accumulation curve for each habitat using all sponges
curve_ATHCO2_calc = specaccum(Calcarea[, 2:94], method = "rarefaction")
curve_ATHCO2_hom = specaccum(Homoscleromorpha[, 2:94], method = "rarefaction") #doesn't work bc there are not any...
curve_ATHCO2_ker = specaccum(Keratosa[, 2:94], method = "rarefaction")
curve_ATHCO2_tetr = specaccum(Tetractinellida[, 2:94], method = "rarefaction")
curve_ATHCO2_sub = specaccum(Suberitida[, 2:94], method = "rarefaction")
curve_ATHCO2_poec = specaccum(Poecilosclerida[, 2:94], method = "rarefaction")
curve_ATHCO2_hap = specaccum(Haplosclerida[, 2:94], method = "rarefaction")
curve_ATHCO2_teth = specaccum(Tethyida[, 2:94], method = "rarefaction")
curve_ATHCO2_ver = specaccum(Verongida[, 2:94], method = "rarefaction")
curve_ATHCO2_Dicty = specaccum(Dictyoceratida[, 2:94], method = "rarefaction")
curve_ATHCO2_Dendr = specaccum(Dendroceratida[, 2:94], method = "rarefaction")
#figure out color palet "display.brewer.pal(n = 8, name = 'RdBu')"
#got the error Error in plot.new() : figure margins too large" so I typed the following command "par(mar=c(5,5,3,1))"" # this fixed the margins and shows axis vallues. 
#plot curve_all first by dates
par(mar=c(5,5,3,1))
plot(curve_ATHCO2_calc, xvar=c("individuals"),ylim=c(0,15),xlim=c(0,800),col="#B2182B", lwd=2, ci.lty=0, ci.col="#F4A582",
     main = "Default: Prettier CI")

#then plot the rest
plot(curve_ATHCO2_hom, add = TRUE,xvar=c("individuals"), col="#2166AC", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI") #col is COLOUR setting, so change it to something else if you #error here because homo has "0

plot(curve_ATHCO2_ker, add = TRUE,xvar=c("individuals"), col="#E1BE6A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_tetr, add = TRUE,xvar=c("individuals"), col="#40B0A6", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_sub, add = TRUE,xvar=c("individuals"), col="#E66100", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_poec, add = TRUE,xvar=c("individuals"), col="#5D3A9B", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_hap, add = TRUE,xvar=c("individuals"), col="#1AFF1A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_teth, add = TRUE,xvar=c("individuals"), col="#994F00", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")


legend(500, 1, legend=c("Calcareous", "Homoscleromorpha"),
       col=c("#B2182B", "#2166AC"), lty=1:1, cex=0.8)

#decided to exclude Dictyo Dendro and Verong since counts were minimal.
```


#SUbset data by TMT="HTHCO2"

```{r setup, include=FALSE}
#Subset data by TMT="HTHCO2" 
HTHCO2physeqpor=subset_samples(physeqpor,TMT=="HTHCO2")
#OTHER MS For ARMS vs ReEf paper phylogenetic analysis etc. 
HTHCO2physeqpormergedDate = merge_samples(HTHCO2physeqpor, "Date")
SD= merge_samples(sample_data(HTHCO2physeqpor), "Date")
print(SD[, "Date"])
print(HTHCO2physeqpormergedDate)
sample_names(HTHCO2physeqpormergedDate)
sample_data(HTHCO2physeqpormergedDate)$Date
identical(SD,sample_data(HTHCO2physeqpormergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqpormergedDate)$Date <- factor(sample_names(HTHCO2physeqpormergedDate))

#Getting abundance
poriferaHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqpormergedDate), "matrix")
View(poriferaHTHCO2_abund_mergedDate)
write.csv(poriferaHTHCO2_abund_mergedDate, "poriferaHTHCO2_abund_mergedDate.csv", row.names = F)

#incorporating abundance changing it to presence absence data
poriferaHTHCO2_abund_mergedDate_abs_pre<-poriferaHTHCO2_abund_mergedDate

poriferaHTHCO2_abund_mergedDate_abs_pre[poriferaHTHCO2_abund_mergedDate_abs_pre >0] <-1 #making a presence absence data frame
View(poriferaHTHCO2_abund_mergedDate_abs_pre)
write.csv(poriferaHTHCO2_abund_mergedDate_abs_pre, "poriferaHTHCO2_abund_mergedDate_abs_pre.csv", row.names = F)

#Getting abundance
poriferaHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqpormergedDate), "matrix")
View(poriferaHTHCO2_abund_mergedunit)
write.csv(poriferaHTHCO2_abund_mergedunit, "poriferaHTHCO2_abund_mergedunit.csv", row.names = F)

###subsetting samples by "HTHCO2"" and then merging sample by "Unit"" 
#which will combine the total diversity of the experiment for this treatment and calculate mean values of diversity#####
HTHCO2physeqpor=subset_samples(physeqpor,TMT=="HTHCO2")
#OTHER MS For ARMS vs ReEf paper phylogenetic analysis etc. 
HTHCO2physeqpormergedunit = merge_samples(HTHCO2physeqpor, "UNIT")
SD= merge_samples(sample_data(HTHCO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(HTHCO2physeqpormergedunit)
sample_names(HTHCO2physeqpormergedunit)
sample_data(HTHCO2physeqpormergedunit)$UNIT
identical(SD,sample_data(HTHCO2physeqpormergedunit))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqpormergedunit)$UNIT <- factor(sample_names(HTHCO2physeqpormergedunit))

#Determine diversity within Porifera HTHCO2.
Richness_HTHCO2physeqpormergedunit<-estimate_richness(HTHCO2physeqpormergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqpormergedunit, "Richness_HTHCO2physeqpormergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqpormergedunit)=sample_names(HTHCO2physeqpormergedunit)
Richness_HTHCO2physeqpormergedunit<-dfRowName(Richness_HTHCO2physeqpormergedunit, name ="UNIT")
View(Richness_HTHCO2physeqpormergedunit)
#need to add a column "TMT" to indicate "HTHCO2"
Richness_HTHCO2physeqpormergedunit<-cbind(Richness_HTHCO2physeqpormergedunit, TMT='HTHCO2')
View(Richness_HTHCO2physeqpormergedunit)
#Calculating mean values of Observed, SHannon, Simpson for the HTHCO2. which will be shown in table xx diversity. 
#Observed
Richness_HTHCO2physeqpormergedunit_sum<-summarySE(Richness_HTHCO2physeqpormergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqpormergedunit_sum)
#Shannon
Shannon_HTHCO2physeqpormergedunit_sum<-summarySE(Richness_HTHCO2physeqpormergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqpormergedunit_sum)
#Simpson
Simpson_HTHCO2physeqpormergedunit_sum<-summarySE(Richness_HTHCO2physeqpormergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqpormergedunit_sum)

#Getting abundance
poriferaHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqpormergedunit), "matrix")
View(poriferaHTHCO2_abund_mergedunit)
write.csv(poriferaHTHCO2_abund_mergedunit, "poriferaHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
poriferaHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(poriferaHTHCO2_abund_mergedunit), data.frame(poriferaHTHCO2_abund_mergedunit, row.names=NULL))
View(poriferaHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poriferaHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
poriferaHTHCO2_abund_mergedunit_asdf<-as_data_frame(poriferaHTHCO2_abund_mergedunit) 
View(poriferaHTHCO2_abund_mergedunit_asdf)
#Sum rows
poriferaHTHCO2_abund_mergedunit_asdf_total<-cbind(poriferaHTHCO2_abund_mergedunit_asdf, total =rowSums(poriferaHTHCO2_abund_mergedunit_asdf)) 
View(poriferaHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poriferaHTHCO2_abund_mergedunit_rn_to_col" to "poriferaHTHCO2_abund_mergedunit_rn_to_col_total"
poriferaHTHCO2_abund_mergedunit_asdf_total$Unit<-poriferaHTHCO2_abund_mergedunit_rn_to_col$Unit
View(poriferaHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 94. 
poriferaHTHCO2_abund_mergedunit_asdf_total_clean<-poriferaHTHCO2_abund_mergedunit_asdf_total[,-(1:93)]
View(poriferaHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
poriferaHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(poriferaHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(poriferaHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(poriferaHTHCO2_abund_mergedunit_asdf_total_clean, "poriferaHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poriferaHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poriferaHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poriferaHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Now calculate diversity within the different sponge subgroups within "HTHCO2physeqpormergedunit"" starting with Class Demospongiae
HTHCO2physeqdemmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Class=="Demospongiae")
Richness_HTHCO2physeqdemmergedunit<-estimate_richness(HTHCO2physeqdemmergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTHCO2physeqdemmergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqdemmergedunit, "Richness_HTHCO2physeqdemmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqdemmergedunit)=sample_names(HTHCO2physeqdemmergedunit)
View(Richness_HTHCO2physeqdemmergedunit)
Richness_HTHCO2physeqdemmergedunit<-dfRowName(Richness_HTHCO2physeqdemmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqdemmergedunit)
#Observed #add column with TMT name HTHCO2
Richness_HTHCO2physeqdemmergedunit <- Richness_HTHCO2physeqdemmergedunit %>%
  add_column(TMT = "HTHCO2")
Richness_HTHCO2physeqdemmergedunit
Richness_HTHCO2physeqdemmergedunit_sum<-summarySE(Richness_HTHCO2physeqdemmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqdemmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqdemmergedunit_sum<-summarySE(Richness_HTHCO2physeqdemmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqdemmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqdemmergedunit_sum<-summarySE(Richness_HTHCO2physeqdemmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqdemmergedunit_sum)

#Getting abundance
demHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqdemmergedunit), "matrix")
View(demHTHCO2_abund_mergedunit)
write.csv(demHTHCO2_abund_mergedunit, "demHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
demHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(demHTHCO2_abund_mergedunit), data.frame(demHTHCO2_abund_mergedunit, row.names=NULL))
View(demHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(demHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
demHTHCO2_abund_mergedunit_asdf<-as_data_frame(demHTHCO2_abund_mergedunit) 
View(demHTHCO2_abund_mergedunit_asdf)
#Sum rows
demHTHCO2_abund_mergedunit_asdf_total<-cbind(demHTHCO2_abund_mergedunit_asdf, total =rowSums(demHTHCO2_abund_mergedunit_asdf)) 
View(demHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "demHTHCO2_abund_mergedunit_rn_to_col" to "demHTHCO2_abund_mergedunit_rn_to_col_total"
demHTHCO2_abund_mergedunit_asdf_total$Unit<-demHTHCO2_abund_mergedunit_rn_to_col$Unit
View(demHTHCO2_abund_mergedunit_asdf_total)
names(demHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 67. 
demHTHCO2_abund_mergedunit_asdf_total_clean<-demHTHCO2_abund_mergedunit_asdf_total[,-(1:66)]
View(demHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
demHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(demHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(demHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(demHTHCO2_abund_mergedunit_asdf_total_clean, "demHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
demHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(demHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(demHTHCO2_abund_mergedunit_asdf_total_clean_sum)


#Now calculate diversity within the different sponge subgroups within "HTHCO2physeqpormergedunit"" starting with Class Homoscleromorpha
HTHCO2physeqhomomergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Class=="Homoscleromorpha")
Richness_HTHCO2physeqhomomergedunit<-estimate_richness(HTHCO2physeqhomomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTHCO2physeqhomomergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqhomomergedunit, "Richness_HTHCO2physeqhomomergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqhomomergedunit)=sample_names(HTHCO2physeqhomomergedunit)
View(Richness_HTHCO2physeqhomomergedunit)
Richness_HTHCO2physeqhomomergedunit<-dfRowName(Richness_HTHCO2physeqhomomergedunit, name ="UNIT")
View(Richness_HTHCO2physeqhomomergedunit)
#Observed #add column with TMT name HTHCO2
Richness_HTHCO2physeqhomomergedunit <- Richness_HTHCO2physeqhomomergedunit %>%
  add_column(TMT = "HTHCO2")
Richness_HTHCO2physeqhomomergedunit
Richness_HTHCO2physeqhomomergedunit_sum<-summarySE(Richness_HTHCO2physeqhomomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqhomomergedunit_sum)
#Shannon
Shannon_HTHCO2physeqhomomergedunit_sum<-summarySE(Richness_HTHCO2physeqhomomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqhomomergedunit_sum)
#Simpson
Simpson_HTHCO2physeqhomomergedunit_sum<-summarySE(Richness_HTHCO2physeqhomomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqhomomergedunit_sum)

#Getting abundance
homoHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqhomomergedunit), "matrix")
View(homoHTHCO2_abund_mergedunit)
write.csv(homoHTHCO2_abund_mergedunit, "homoHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
homoHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(homoHTHCO2_abund_mergedunit), data.frame(homoHTHCO2_abund_mergedunit, row.names=NULL))
View(homoHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(homoHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
homoHTHCO2_abund_mergedunit_asdf<-as_data_frame(homoHTHCO2_abund_mergedunit) 
View(homoHTHCO2_abund_mergedunit_asdf)
#Sum rows
homoHTHCO2_abund_mergedunit_asdf_total<-cbind(homoHTHCO2_abund_mergedunit_asdf, total =rowSums(homoHTHCO2_abund_mergedunit_asdf)) 
View(homoHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "homoHTHCO2_abund_mergedunit_rn_to_col" to "homoHTHCO2_abund_mergedunit_rn_to_col_total"
homoHTHCO2_abund_mergedunit_asdf_total$Unit<-homoHTHCO2_abund_mergedunit_rn_to_col$Unit
View(homoHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
homoHTHCO2_abund_mergedunit_asdf_total_clean<-homoHTHCO2_abund_mergedunit_asdf_total[,-(1:8)]
View(homoHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
homoHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(homoHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(homoHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(homoHTHCO2_abund_mergedunit_asdf_total_clean, "homoHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
homoHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(homoHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(homoHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Dictyoceratida
#Now calculate diversity within the different sponge subgroups within "HTHCO2physeqpormergedunit"" starting with Order Dictyoceratida
HTHCO2physeqdictyomergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Dictyoceratida")
Richness_HTHCO2physeqdictyomergedunit<-estimate_richness(HTHCO2physeqdictyomergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTHCO2physeqdictyomergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqdictyomergedunit, "Richness_HTHCO2physeqdictyomergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqdictyomergedunit)=sample_names(HTHCO2physeqdictyomergedunit)
View(Richness_HTHCO2physeqdictyomergedunit)
Richness_HTHCO2physeqdictyomergedunit<-dfRowName(Richness_HTHCO2physeqdictyomergedunit, name ="UNIT")
View(Richness_HTHCO2physeqdictyomergedunit)
#Observed #add column with TMT name HTHCO2
Richness_HTHCO2physeqdictyomergedunit <- Richness_HTHCO2physeqdictyomergedunit %>%
  add_column(TMT = "HTHCO2")
Richness_HTHCO2physeqdictyomergedunit
Richness_HTHCO2physeqdictyomergedunit_sum<-summarySE(Richness_HTHCO2physeqdictyomergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqdictyomergedunit_sum)
#Shannon
Shannon_HTHCO2physeqdictyomergedunit_sum<-summarySE(Richness_HTHCO2physeqdictyomergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqdictyomergedunit_sum)
#Simpson
Simpson_HTHCO2physeqdictyomergedunit_sum<-summarySE(Richness_HTHCO2physeqdictyomergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqdictyomergedunit_sum)
#Getting abundance
dictyoHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqdictyomergedunit), "matrix")
View(dictyoHTHCO2_abund_mergedunit)
write.csv(dictyoHTHCO2_abund_mergedunit, "dictyoHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dictyoHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(dictyoHTHCO2_abund_mergedunit), data.frame(dictyoHTHCO2_abund_mergedunit, row.names=NULL))
View(dictyoHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dictyoHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dictyoHTHCO2_abund_mergedunit_asdf<-as_data_frame(dictyoHTHCO2_abund_mergedunit) 
View(dictyoHTHCO2_abund_mergedunit_asdf)
#Sum rows
dictyoHTHCO2_abund_mergedunit_asdf_total<-cbind(dictyoHTHCO2_abund_mergedunit_asdf, total =rowSums(dictyoHTHCO2_abund_mergedunit_asdf)) 
View(dictyoHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dictyoHTHCO2_abund_mergedunit_rn_to_col" to "dictyoHTHCO2_abund_mergedunit_rn_to_col_total"
dictyoHTHCO2_abund_mergedunit_asdf_total$Unit<-dictyoHTHCO2_abund_mergedunit_rn_to_col$Unit
View(dictyoHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dictyoHTHCO2_abund_mergedunit_asdf_total_clean<-dictyoHTHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dictyoHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
dictyoHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(dictyoHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(dictyoHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(dictyoHTHCO2_abund_mergedunit_asdf_total_clean, "dictyoHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dictyoHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dictyoHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dictyoHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Dendroceratida
#Now calculate diversity within the different sponge subgroups within "HTHCO2physeqpormergedunit"" starting with Order Dendroceratida
HTHCO2physeqdendromergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Dendroceratida")
Richness_HTHCO2physeqdendromergedunit<-estimate_richness(HTHCO2physeqdendromergedunit, measures = c("Observed", "Shannon", "Simpson"))
View(Richness_HTHCO2physeqdendromergedunit)
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqdendromergedunit, "Richness_HTHCO2physeqdendromergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqdendromergedunit)=sample_names(HTHCO2physeqdendromergedunit)
View(Richness_HTHCO2physeqdendromergedunit)
Richness_HTHCO2physeqdendromergedunit<-dfRowName(Richness_HTHCO2physeqdendromergedunit, name ="UNIT")
View(Richness_HTHCO2physeqdendromergedunit)
#Observed #add column with TMT name HTHCO2
Richness_HTHCO2physeqdendromergedunit <- Richness_HTHCO2physeqdendromergedunit %>%
  add_column(TMT = "HTHCO2")
Richness_HTHCO2physeqdendromergedunit_sum<-summarySE(Richness_HTHCO2physeqdendromergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqdendromergedunit_sum)
#Shannon
Shannon_HTHCO2physeqdendromergedunit_sum<-summarySE(Richness_HTHCO2physeqdendromergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqdendromergedunit_sum)
#Simpson
Simpson_HTHCO2physeqdendromergedunit_sum<-summarySE(Richness_HTHCO2physeqdendromergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqdendromergedunit_sum)
#Getting abundance
dendroHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqdendromergedunit), "matrix")
View(dendroHTHCO2_abund_mergedunit)
write.csv(dendroHTHCO2_abund_mergedunit, "dendroHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
dendroHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(dendroHTHCO2_abund_mergedunit), data.frame(dendroHTHCO2_abund_mergedunit, row.names=NULL))
View(dendroHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(dendroHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
dendroHTHCO2_abund_mergedunit_asdf<-as_data_frame(dendroHTHCO2_abund_mergedunit) 
View(dendroHTHCO2_abund_mergedunit_asdf)
#Sum rows
dendroHTHCO2_abund_mergedunit_asdf_total<-cbind(dendroHTHCO2_abund_mergedunit_asdf, total =rowSums(dendroHTHCO2_abund_mergedunit_asdf)) 
View(dendroHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "dendroHTHCO2_abund_mergedunit_rn_to_col" to "dendroHTHCO2_abund_mergedunit_rn_to_col_total"
dendroHTHCO2_abund_mergedunit_asdf_total$Unit<-dendroHTHCO2_abund_mergedunit_rn_to_col$Unit
View(dendroHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
dendroHTHCO2_abund_mergedunit_asdf_total_clean<-dendroHTHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(dendroHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
dendroHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(dendroHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(dendroHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(dendroHTHCO2_abund_mergedunit_asdf_total_clean, "dendroHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
dendroHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(dendroHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(dendroHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Calcarea
HTHCO2physeqcalcmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Class=="Calcarea")
Richness_HTHCO2physeqcalcmergedunit<-estimate_richness(HTHCO2physeqcalcmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqcalcmergedunit, "Richness_HTHCO2physeqcalcmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqcalcmergedunit)=sample_names(HTHCO2physeqcalcmergedunit)
Richness_HTHCO2physeqcalcmergedunit<-dfRowName(Richness_HTHCO2physeqcalcmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqcalcmergedunit)
Richness_HTHCO2physeqcalcmergedunit<-cbind(Richness_HTHCO2physeqcalcmergedunit, TMT='HTHCO2')
View(Richness_HTHCO2physeqcalcmergedunit)
#Observed
Richness_HTHCO2physeqcalcmergedunit_sum<-summarySE(Richness_HTHCO2physeqcalcmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqcalcmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqcalcmergedunit_sum<-summarySE(Richness_HTHCO2physeqcalcmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqcalcmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqcalcmergedunit_sum<-summarySE(Richness_HTHCO2physeqcalcmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqcalcmergedunit_sum)
#Getting abundance
calcHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqcalcmergedunit), "matrix")
View(calcHTHCO2_abund_mergedunit)
write.csv(calcHTHCO2_abund_mergedunit, "calcHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
calcHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(calcHTHCO2_abund_mergedunit), data.frame(calcHTHCO2_abund_mergedunit, row.names=NULL))
View(calcHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(calcHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
calcHTHCO2_abund_mergedunit_asdf<-as_data_frame(calcHTHCO2_abund_mergedunit) 
View(calcHTHCO2_abund_mergedunit_asdf)
#Sum rows
calcHTHCO2_abund_mergedunit_asdf_total<-cbind(calcHTHCO2_abund_mergedunit_asdf, total =rowSums(calcHTHCO2_abund_mergedunit_asdf)) 
View(calcHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "calcHTHCO2_abund_mergedunit_rn_to_col" to "calcHTHCO2_abund_mergedunit_rn_to_col_total"
calcHTHCO2_abund_mergedunit_asdf_total$Unit<-calcHTHCO2_abund_mergedunit_rn_to_col$Unit
View(calcHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
calcHTHCO2_abund_mergedunit_asdf_total_clean<-calcHTHCO2_abund_mergedunit_asdf_total[,-(1:19)]
View(calcHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
calcHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(calcHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(calcHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(calcHTHCO2_abund_mergedunit_asdf_total_clean, "calcHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
calcHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(calcHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(calcHTHCO2_abund_mergedunit_asdf_total_clean_sum)


#Keratosa
HTHCO2physeqkermergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Subclass=="Keratosa")
Richness_HTHCO2physeqkermergedunit<-estimate_richness(HTHCO2physeqkermergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqkermergedunit, "Richness_HTHCO2physeqkermergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqkermergedunit)=sample_names(HTHCO2physeqkermergedunit)
Richness_HTHCO2physeqkermergedunit<-dfRowName(Richness_HTHCO2physeqkermergedunit, name ="UNIT")
View(Richness_HTHCO2physeqkermergedunit)
Richness_HTHCO2physeqkermergedunit<-cbind(Richness_HTHCO2physeqkermergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqkermergedunit_sum<-summarySE(Richness_HTHCO2physeqkermergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqkermergedunit_sum)
#Shannon
Shannon_HTHCO2physeqkermergedunit_sum<-summarySE(Richness_HTHCO2physeqkermergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqkermergedunit_sum)
#Simpson
Simpson_HTHCO2physeqkermergedunit_sum<-summarySE(Richness_HTHCO2physeqkermergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqkermergedunit_sum)
#Getting abundance
kerHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqkermergedunit), "matrix")
View(kerHTHCO2_abund_mergedunit)
write.csv(kerHTHCO2_abund_mergedunit, "kerHTHCO2_abund_mergedunit.csv", row.names = F)
#Convert Rowname into first column 
kerHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(kerHTHCO2_abund_mergedunit), data.frame(kerHTHCO2_abund_mergedunit, row.names=NULL))
View(kerHTHCO2_abund_mergedunit)
#Change rowname to "Unit"
colnames(kerHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#convert to dataframe
kerHTHCO2_abund_mergedunit_asdf<-as_data_frame(kerHTHCO2_abund_mergedunit) 
View(kerHTHCO2_abund_mergedunit_asdf)
#Sum rows
kerHTHCO2_abund_mergedunit_asdf_total<-cbind(kerHTHCO2_abund_mergedunit_asdf, total =rowSums(kerHTHCO2_abund_mergedunit_asdf)) 
View(kerHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "kerHTHCO2_abund_mergedunit_rn_to_col" to "kerHTHCO2_abund_mergedunit_rn_to_col_total"
kerHTHCO2_abund_mergedunit_asdf_total$Unit<-kerHTHCO2_abund_mergedunit_rn_to_col$Unit
View(kerHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 93. 
kerHTHCO2_abund_mergedunit_asdf_total_clean<-kerHTHCO2_abund_mergedunit_asdf_total[,-(1:6)]
View(kerHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
kerHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(kerHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(kerHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(kerHTHCO2_abund_mergedunit_asdf_total_clean, "kerHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
kerHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(kerHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(kerHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Tetractinellida
HTHCO2physeqtetrmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Tetractinellida")
Richness_HTHCO2physeqtetrmergedunit<-estimate_richness(HTHCO2physeqtetrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqtetrmergedunit, "Richness_HTHCO2physeqtetrmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqtetrmergedunit)=sample_names(HTHCO2physeqtetrmergedunit)
Richness_HTHCO2physeqtetrmergedunit<-dfRowName(Richness_HTHCO2physeqtetrmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqtetrmergedunit)
Richness_HTHCO2physeqtetrmergedunit<-cbind(Richness_HTHCO2physeqtetrmergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqtetrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtetrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqtetrmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqtetrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtetrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqtetrmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqtetrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtetrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqtetrmergedunit_sum)
#Getting abundance
tetrHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqtetrmergedunit), "matrix")
View(tetrHTHCO2_abund_mergedunit)
write.csv(tetrHTHCO2_abund_mergedunit, "tetrHTHCO2_abund_mergedunit.csv", row.names = F)
#Contetrt Rowname into first column 
tetrHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(tetrHTHCO2_abund_mergedunit), data.frame(tetrHTHCO2_abund_mergedunit, row.names=NULL))
View(tetrHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tetrHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contetrt to dataframe
tetrHTHCO2_abund_mergedunit_asdf<-as_data_frame(tetrHTHCO2_abund_mergedunit) 
View(tetrHTHCO2_abund_mergedunit_asdf)
#Sum rows
tetrHTHCO2_abund_mergedunit_asdf_total<-cbind(tetrHTHCO2_abund_mergedunit_asdf, total =rowSums(tetrHTHCO2_abund_mergedunit_asdf)) 
View(tetrHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tetrHTHCO2_abund_mergedunit_rn_to_col" to "tetrHTHCO2_abund_mergedunit_rn_to_col_total"
tetrHTHCO2_abund_mergedunit_asdf_total$Unit<-tetrHTHCO2_abund_mergedunit_rn_to_col$Unit
View(tetrHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 3. 
tetrHTHCO2_abund_mergedunit_asdf_total_clean<-tetrHTHCO2_abund_mergedunit_asdf_total[,-(1:3)]
View(tetrHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
tetrHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(tetrHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(tetrHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(tetrHTHCO2_abund_mergedunit_asdf_total_clean, "tetrHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tetrHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tetrHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tetrHTHCO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Poecilosclerida
HTHCO2physeqpoecrmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Poecilosclerida")
Richness_HTHCO2physeqpoecrmergedunit<-estimate_richness(HTHCO2physeqpoecrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqpoecrmergedunit, "Richness_HTHCO2physeqpoecrmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqpoecrmergedunit)=sample_names(HTHCO2physeqpoecrmergedunit)
Richness_HTHCO2physeqpoecrmergedunit<-dfRowName(Richness_HTHCO2physeqpoecrmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqpoecrmergedunit)
Richness_HTHCO2physeqpoecrmergedunit<-cbind(Richness_HTHCO2physeqpoecrmergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTHCO2physeqpoecrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqpoecrmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTHCO2physeqpoecrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqpoecrmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqpoecrmergedunit_sum<-summarySE(Richness_HTHCO2physeqpoecrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqpoecrmergedunit_sum)
#Getting abundance
poecrHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqpoecrmergedunit), "matrix")
View(poecrHTHCO2_abund_mergedunit)
write.csv(poecrHTHCO2_abund_mergedunit, "poecrHTHCO2_abund_mergedunit.csv", row.names = F)
#Conpoecrt Rowname into first column 
poecrHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(poecrHTHCO2_abund_mergedunit), data.frame(poecrHTHCO2_abund_mergedunit, row.names=NULL))
View(poecrHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(poecrHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conpoecrt to dataframe
poecrHTHCO2_abund_mergedunit_asdf<-as_data_frame(poecrHTHCO2_abund_mergedunit) 
View(poecrHTHCO2_abund_mergedunit_asdf)
#Sum rows
poecrHTHCO2_abund_mergedunit_asdf_total<-cbind(poecrHTHCO2_abund_mergedunit_asdf, total =rowSums(poecrHTHCO2_abund_mergedunit_asdf)) 
View(poecrHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "poecrHTHCO2_abund_mergedunit_rn_to_col" to "poecrHTHCO2_abund_mergedunit_rn_to_col_total"
poecrHTHCO2_abund_mergedunit_asdf_total$Unit<-poecrHTHCO2_abund_mergedunit_rn_to_col$Unit
View(poecrHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
poecrHTHCO2_abund_mergedunit_asdf_total_clean<-poecrHTHCO2_abund_mergedunit_asdf_total[,-(1:12)]
View(poecrHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
poecrHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(poecrHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(poecrHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(poecrHTHCO2_abund_mergedunit_asdf_total_clean, "poecrHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
poecrHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(poecrHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(poecrHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Tethyida
HTHCO2physeqtethrmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Tethyida")
Richness_HTHCO2physeqtethrmergedunit<-estimate_richness(HTHCO2physeqtethrmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqtethrmergedunit, "Richness_HTHCO2physeqtethrmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqtethrmergedunit)=sample_names(HTHCO2physeqtethrmergedunit)
Richness_HTHCO2physeqtethrmergedunit<-dfRowName(Richness_HTHCO2physeqtethrmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqtethrmergedunit)
Richness_HTHCO2physeqtethrmergedunit<-cbind(Richness_HTHCO2physeqtethrmergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqtethrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtethrmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqtethrmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqtethrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtethrmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqtethrmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqtethrmergedunit_sum<-summarySE(Richness_HTHCO2physeqtethrmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqtethrmergedunit_sum)
#Getting abundance
tethrHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqtethrmergedunit), "matrix")
View(tethrHTHCO2_abund_mergedunit)
write.csv(tethrHTHCO2_abund_mergedunit, "tethrHTHCO2_abund_mergedunit.csv", row.names = F)
#Contethrt Rowname into first column 
tethrHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(tethrHTHCO2_abund_mergedunit), data.frame(tethrHTHCO2_abund_mergedunit, row.names=NULL))
View(tethrHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(tethrHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#contethrt to dataframe
tethrHTHCO2_abund_mergedunit_asdf<-as_data_frame(tethrHTHCO2_abund_mergedunit) 
View(tethrHTHCO2_abund_mergedunit_asdf)
#Sum rows
tethrHTHCO2_abund_mergedunit_asdf_total<-cbind(tethrHTHCO2_abund_mergedunit_asdf, total =rowSums(tethrHTHCO2_abund_mergedunit_asdf)) 
View(tethrHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "tethrHTHCO2_abund_mergedunit_rn_to_col" to "tethrHTHCO2_abund_mergedunit_rn_to_col_total"
tethrHTHCO2_abund_mergedunit_asdf_total$Unit<-tethrHTHCO2_abund_mergedunit_rn_to_col$Unit
View(tethrHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 12. 
tethrHTHCO2_abund_mergedunit_asdf_total_clean<-tethrHTHCO2_abund_mergedunit_asdf_total[,-(1:6)]
View(tethrHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
tethrHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(tethrHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(tethrHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(tethrHTHCO2_abund_mergedunit_asdf_total_clean, "tethrHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
tethrHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(tethrHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(tethrHTHCO2_abund_mergedunit_asdf_total_clean_sum)


#Order==Haplosclerida
HTHCO2physeqhaprmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Haplosclerida")
Richness_HTHCO2physeqhaprmergedunit<-estimate_richness(HTHCO2physeqhaprmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqhaprmergedunit, "Richness_HTHCO2physeqhaprmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqhaprmergedunit)=sample_names(HTHCO2physeqhaprmergedunit)
Richness_HTHCO2physeqhaprmergedunit<-dfRowName(Richness_HTHCO2physeqhaprmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqhaprmergedunit)
Richness_HTHCO2physeqhaprmergedunit<-cbind(Richness_HTHCO2physeqhaprmergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqhaprmergedunit_sum<-summarySE(Richness_HTHCO2physeqhaprmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqhaprmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqhaprmergedunit_sum<-summarySE(Richness_HTHCO2physeqhaprmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqhaprmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqhaprmergedunit_sum<-summarySE(Richness_HTHCO2physeqhaprmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqhaprmergedunit_sum)
#Getting abundance
haprHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqhaprmergedunit), "matrix")
View(haprHTHCO2_abund_mergedunit)
write.csv(haprHTHCO2_abund_mergedunit, "haprHTHCO2_abund_mergedunit.csv", row.names = F)
#Conhaprt Rowname into first column 
haprHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(haprHTHCO2_abund_mergedunit), data.frame(haprHTHCO2_abund_mergedunit, row.names=NULL))
View(haprHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(haprHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#conhaprt to dataframe
haprHTHCO2_abund_mergedunit_asdf<-as_data_frame(haprHTHCO2_abund_mergedunit) 
View(haprHTHCO2_abund_mergedunit_asdf)
#Sum rows
haprHTHCO2_abund_mergedunit_asdf_total<-cbind(haprHTHCO2_abund_mergedunit_asdf, total =rowSums(haprHTHCO2_abund_mergedunit_asdf)) 
View(haprHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "haprHTHCO2_abund_mergedunit_rn_to_col" to "haprHTHCO2_abund_mergedunit_rn_to_col_total"
haprHTHCO2_abund_mergedunit_asdf_total$Unit<-haprHTHCO2_abund_mergedunit_rn_to_col$Unit
View(haprHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 23. 
haprHTHCO2_abund_mergedunit_asdf_total_clean<-haprHTHCO2_abund_mergedunit_asdf_total[,-(1:22)]
View(haprHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
haprHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(haprHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(haprHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(haprHTHCO2_abund_mergedunit_asdf_total_clean, "haprHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
haprHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(haprHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(haprHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Order==Suberitida
HTHCO2physeqsubmergedunit<-subset_taxa(HTHCO2physeqpormergedunit, Order=="Suberitida")
Richness_HTHCO2physeqsubmergedunit<-estimate_richness(HTHCO2physeqsubmergedunit, measures = c("Observed", "Shannon", "Simpson"))
#Export richness values for t test analysis. 
write.csv(Richness_HTHCO2physeqsubmergedunit, "Richness_HTHCO2physeqsubmergedunit.csv", row.names = F)
row.names(Richness_HTHCO2physeqsubmergedunit)=sample_names(HTHCO2physeqsubmergedunit)
Richness_HTHCO2physeqsubmergedunit<-dfRowName(Richness_HTHCO2physeqsubmergedunit, name ="UNIT")
View(Richness_HTHCO2physeqsubmergedunit)
Richness_HTHCO2physeqsubmergedunit<-cbind(Richness_HTHCO2physeqsubmergedunit, TMT='HTHCO2')
#Observed
Richness_HTHCO2physeqsubmergedunit_sum<-summarySE(Richness_HTHCO2physeqsubmergedunit, measurevar = "Observed", groupvars=c("TMT"))
View(Richness_HTHCO2physeqsubmergedunit_sum)
#Shannon
Shannon_HTHCO2physeqsubmergedunit_sum<-summarySE(Richness_HTHCO2physeqsubmergedunit, measurevar = "Shannon", groupvars=c("TMT"))
View(Shannon_HTHCO2physeqsubmergedunit_sum)
#Simpson
Simpson_HTHCO2physeqsubmergedunit_sum<-summarySE(Richness_HTHCO2physeqsubmergedunit, measurevar = "Simpson", groupvars=c("TMT"))
View(Simpson_HTHCO2physeqsubmergedunit_sum)
#Getting abundance
subHTHCO2_abund_mergedunit=as(otu_table(HTHCO2physeqsubmergedunit), "matrix")
View(subHTHCO2_abund_mergedunit)
write.csv(subHTHCO2_abund_mergedunit, "subHTHCO2_abund_mergedunit.csv", row.names = F)
#Consubt Rowname into first column 
subHTHCO2_abund_mergedunit_rn_to_col <- cbind(rownames(subHTHCO2_abund_mergedunit), data.frame(subHTHCO2_abund_mergedunit, row.names=NULL))
View(subHTHCO2_abund_mergedunit_rn_to_col)
#Change rowname to "Unit"
colnames(subHTHCO2_abund_mergedunit_rn_to_col)[1] <- "Unit"
#consubt to dataframe
subHTHCO2_abund_mergedunit_asdf<-as_data_frame(subHTHCO2_abund_mergedunit) 
View(subHTHCO2_abund_mergedunit_asdf)
#Sum rows
subHTHCO2_abund_mergedunit_asdf_total<-cbind(subHTHCO2_abund_mergedunit_asdf, total =rowSums(subHTHCO2_abund_mergedunit_asdf)) 
View(subHTHCO2_abund_mergedunit_asdf_total)
#add column "Unit" from "subHTHCO2_abund_mergedunit_rn_to_col" to "subHTHCO2_abund_mergedunit_rn_to_col_total"
subHTHCO2_abund_mergedunit_asdf_total$Unit<-subHTHCO2_abund_mergedunit_rn_to_col$Unit
View(subHTHCO2_abund_mergedunit_asdf_total)
#Eliminate column 1 through 23. 
subHTHCO2_abund_mergedunit_asdf_total_clean<-subHTHCO2_abund_mergedunit_asdf_total[,-(1:15)]
View(subHTHCO2_abund_mergedunit_asdf_total_clean)
#need to add a column "TMT" to indicate "HTHCO2"
subHTHCO2_abund_mergedunit_asdf_total_clean<-cbind(subHTHCO2_abund_mergedunit_asdf_total_clean, TMT='HTHCO2')
View(subHTHCO2_abund_mergedunit_asdf_total_clean)
write.csv(subHTHCO2_abund_mergedunit_asdf_total_clean, "subHTHCO2_abund_mergedunit_asdf_total_clean.csv", row.names = F)
#SUmmarySE
subHTHCO2_abund_mergedunit_asdf_total_clean_sum<-summarySE(subHTHCO2_abund_mergedunit_asdf_total_clean, measurevar = "total", groupvars=c("TMT"))
View(subHTHCO2_abund_mergedunit_asdf_total_clean_sum)

#Rarefraction curve preliminary steps:Subset data within "HTHCO2DATE"to obtain the different subsponge groups.
#remember that when subsetting taxonomic assignments you must use "subset_taxa" starting with 
#Homoscleromorpha
HTHCO2physeqhomo=subset_taxa(HTHCO2physeqpor,Class=="Homoscleromorpha")
HTHCO2physeqhomomergedDate = merge_samples(HTHCO2physeqhomo, "Date")
SD= merge_samples(sample_data(HTHCO2physeqhomo), "Date")
print(SD[, "Date"])
print(HTHCO2physeqhomomergedDate)
sample_names(HTHCO2physeqhomomergedDate)
sample_data(HTHCO2physeqhomomergedDate)$Date
identical(SD,sample_data(HTHCO2physeqhomomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqhomomergedDate)$Date <- factor(sample_names(HTHCO2physeqhomomergedDate))
#Getting abundance
homoHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqhomomergedDate), "matrix")
View(homoHTHCO2_abund_mergedDate)
#adding column at the end that shows TMT=HTHCO2
homoHTHCO2_abund_mergedDate<-cbind(homoHTHCO2_abund_mergedDate, TMT='HTHCO2')
head(homoHTHCO2_abund_mergedDate)
write.csv(homoHTHCO2_abund_mergedDate, "homoHTHCO2_abund_mergedDate.csv", row.names = F)

#Demospongiae
HTHCO2physeqdem=subset_taxa(HTHCO2physeqpor,Class=="Demospongiae")
HTHCO2physeqdemmergedDate = merge_samples(HTHCO2physeqdem, "Date")
SD= merge_samples(sample_data(HTHCO2physeqdem), "Date")
print(SD[, "Date"])
print(HTHCO2physeqdemmergedDate)
sample_names(HTHCO2physeqdemmergedDate)
sample_data(HTHCO2physeqdemmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqdemmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqdemmergedDate)$Date <- factor(sample_names(HTHCO2physeqdemmergedDate))
#Getting abundance
demHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqdemmergedDate), "matrix")
View(demHTHCO2_abund_mergedDate)
#adding column at the end that shows TMT=HTHCO2
demHTHCO2_abund_mergedDate<-cbind(demHTHCO2_abund_mergedDate, TMT='HTHCO2')
write.csv(demHTHCO2_abund_mergedDate, "demHTHCO2_abund_mergedDate.csv", row.names = F)


#Dictyoceratida
HTHCO2physeqdictyo=subset_taxa(HTHCO2physeqpor,Order=="Dictyoceratida")
HTHCO2physeqdictyomergedDate = merge_samples(HTHCO2physeqdictyo, "Date")
SD= merge_samples(sample_data(HTHCO2physeqdictyo), "Date")
print(SD[, "Date"])
print(HTHCO2physeqdictyomergedDate)
sample_names(HTHCO2physeqdictyomergedDate)
sample_data(HTHCO2physeqdictyomergedDate)$Date
identical(SD,sample_data(HTHCO2physeqdictyomergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqdictyomergedDate)$Date <- factor(sample_names(HTHCO2physeqdictyomergedDate))
#Getting abundance
dictyoHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqdictyomergedDate), "matrix")
View(dictyoHTHCO2_abund_mergedDate)
write.csv(dictyoHTHCO2_abund_mergedDate, "dictyoHTHCO2_abund_mergedDate.csv", row.names = F)

#Dendroceratida
HTHCO2physeqdendro=subset_taxa(HTHCO2physeqpor,Order=="Dendroceratida")
HTHCO2physeqdendromergedDate = merge_samples(HTHCO2physeqdendro, "Date")
SD= merge_samples(sample_data(HTHCO2physeqdendro), "Date")
print(SD[, "Date"])
print(HTHCO2physeqdendromergedDate)
sample_names(HTHCO2physeqdendromergedDate)
sample_data(HTHCO2physeqdendromergedDate)$Date
identical(SD,sample_data(HTHCO2physeqdendromergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqdendromergedDate)$Date <- factor(sample_names(HTHCO2physeqdendromergedDate))
#Getting abundance
dendroHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqdendromergedDate), "matrix")
View(dendroHTHCO2_abund_mergedDate)
write.csv(dendroHTHCO2_abund_mergedDate, "dendroHTHCO2_abund_mergedDate.csv", row.names = F)


###Subset with Clathrinida 
HTHCO2physeqclat=subset_taxa(HTHCO2physeqpor,Order=="Clathrinida")
HTHCO2physeqclatrmergedDate = merge_samples(HTHCO2physeqclat, "Date")
SD= merge_samples(sample_data(HTHCO2physeqclat), "Date")
print(SD[, "Date"])
print(HTHCO2physeqclatrmergedDate)
sample_names(HTHCO2physeqclatrmergedDate)
sample_data(HTHCO2physeqclatrmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqclatrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqclatrmergedDate)$Date <- factor(sample_names(HTHCO2physeqclatrmergedDate))
#Getting abundance
clatHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqclatrmergedDate), "matrix")
View(clatHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
clatHTHCO2_abund_mergedDate_abs_pre<-clatHTHCO2_abund_mergedDate
write.csv(clatHTHCO2_abund_mergedDate, "clatHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Leucosolenida
HTHCO2physeqleuc=subset_taxa(HTHCO2physeqpor,Order=="Leucosolenida")
HTHCO2physeqleucmergedDate = merge_samples(HTHCO2physeqleuc, "Date")
SD= merge_samples(sample_data(HTHCO2physeqleuc), "Date")
print(SD[, "Date"])
print(HTHCO2physeqleucmergedDate)
sample_names(HTHCO2physeqleucmergedDate)
sample_data(HTHCO2physeqleucmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqleucmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqleucmergedDate)$Date <- factor(sample_names(HTHCO2physeqleucmergedDate))
#Getting abundance
leucHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqleucmergedDate), "matrix")
View(leucHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
leucHTHCO2_abund_mergedDate_abs_pre<-leucHTHCO2_abund_mergedDate
write.csv(leucHTHCO2_abund_mergedDate, "leucHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Class: Calcarea
HTHCO2physeqcalc=subset_taxa(HTHCO2physeqpor,Class=="Calcarea")
HTHCO2physeqcalcmergedDate = merge_samples(HTHCO2physeqcalc, "Date")
SD= merge_samples(sample_data(HTHCO2physeqcalc), "Date")
print(SD[, "Date"])
print(HTHCO2physeqcalcmergedDate)
sample_names(HTHCO2physeqcalcmergedDate)
sample_data(HTHCO2physeqcalcmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqcalcmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqcalcmergedDate)$Date <- factor(sample_names(HTHCO2physeqcalcmergedDate))
#Getting abundance
calcHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqcalcmergedDate), "matrix")
View(calcHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
calcHTHCO2_abund_mergedDate_abs_pre<-calcHTHCO2_abund_mergedDate
#adding column at the end that shows class=Demospongiae
calcHTHCO2_abund_mergedDate<-cbind(calcHTHCO2_abund_mergedDate, TMT='HTHCO2')
write.csv(calcHTHCO2_abund_mergedDate, "calcHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Keratosa 
HTHCO2physeqker=subset_taxa(HTHCO2physeqpor,Subclass=="Keratosa")
HTHCO2physeqkermergedDate = merge_samples(HTHCO2physeqker, "Date")
SD= merge_samples(sample_data(HTHCO2physeqker), "Date")
print(SD[, "Date"])
print(HTHCO2physeqkermergedDate)
sample_names(HTHCO2physeqkermergedDate)
sample_data(HTHCO2physeqkermergedDate)$Date
identical(SD,sample_data(HTHCO2physeqkermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqkermergedDate)$Date <- factor(sample_names(HTHCO2physeqkermergedDate))
#Getting abundance
kerHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqkermergedDate), "matrix")
View(kerHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
kerHTHCO2_abund_mergedDate_abs_pre<-kerHTHCO2_abund_mergedDate
write.csv(kerHTHCO2_abund_mergedDate, "kerHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Tetractinellida 
HTHCO2physeqtetr=subset_taxa(HTHCO2physeqpor,Order=="Tetractinellida")
HTHCO2physeqtetrmergedDate = merge_samples(HTHCO2physeqtetr, "Date")
SD= merge_samples(sample_data(HTHCO2physeqtetr), "Date")
print(SD[, "Date"])
print(HTHCO2physeqtetrmergedDate)
sample_names(HTHCO2physeqtetrmergedDate)
sample_data(HTHCO2physeqtetrmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqtetrmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqtetrmergedDate)$Date <- factor(sample_names(HTHCO2physeqtetrmergedDate))
#Getting abundance
tetrHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqtetrmergedDate), "matrix")
View(tetrHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tetrHTHCO2_abund_mergedDate_abs_pre<-tetrHTHCO2_abund_mergedDate
write.csv(tetrHTHCO2_abund_mergedDate, "tetrHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Suberitida 
HTHCO2physeqsub=subset_taxa(HTHCO2physeqpor,Order=="Suberitida")
HTHCO2physeqsubmergedDate = merge_samples(HTHCO2physeqsub, "Date")
SD= merge_samples(sample_data(HTHCO2physeqsub), "Date")
print(SD[, "Date"])
print(HTHCO2physeqsubmergedDate)
sample_names(HTHCO2physeqsubmergedDate)
sample_data(HTHCO2physeqsubmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqsubmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqsubmergedDate)$Date <- factor(sample_names(HTHCO2physeqsubmergedDate))
#Getting abundance
subHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqsubmergedDate), "matrix")
View(subHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
subHTHCO2_abund_mergedDate_abs_pre<-subHTHCO2_abund_mergedDate
write.csv(subHTHCO2_abund_mergedDate, "subHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Poecilosclerida
HTHCO2physeqpoec=subset_taxa(HTHCO2physeqpor,Order=="Poecilosclerida")
HTHCO2physeqpoecmergedDate = merge_samples(HTHCO2physeqpoec, "Date")
SD= merge_samples(sample_data(HTHCO2physeqpoec), "Date")
print(SD[, "Date"])
print(HTHCO2physeqpoecmergedDate)
sample_names(HTHCO2physeqpoecmergedDate)
sample_data(HTHCO2physeqpoecmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqpoecmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqpoecmergedDate)$Date <- factor(sample_names(HTHCO2physeqpoecmergedDate))
#Getting abundance
poecHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqpoecmergedDate), "matrix")
View(poecHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
poecHTHCO2_abund_mergedDate_abs_pre<-poecHTHCO2_abund_mergedDate
write.csv(poecHTHCO2_abund_mergedDate, "poecHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Haplosclerida
HTHCO2physeqhap=subset_taxa(HTHCO2physeqpor,Order=="Haplosclerida")
HTHCO2physeqhapmergedDate = merge_samples(HTHCO2physeqhap, "Date")
SD= merge_samples(sample_data(HTHCO2physeqhap), "Date")
print(SD[, "Date"])
print(HTHCO2physeqhapmergedDate)
sample_names(HTHCO2physeqhapmergedDate)
sample_data(HTHCO2physeqhapmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqhapmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqhapmergedDate)$Date <- factor(sample_names(HTHCO2physeqhapmergedDate))
#Getting abundance
hapHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqhapmergedDate), "matrix")
View(hapHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
hapHTHCO2_abund_mergedDate_abs_pre<-hapHTHCO2_abund_mergedDate
write.csv(hapHTHCO2_abund_mergedDate, "hapHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Tethyida
HTHCO2physeqteth=subset_taxa(HTHCO2physeqpor,Order=="Tethyida")
HTHCO2physeqtethmergedDate = merge_samples(HTHCO2physeqteth, "Date")
SD= merge_samples(sample_data(HTHCO2physeqteth), "Date")
print(SD[, "Date"])
print(HTHCO2physeqtethmergedDate)
sample_names(HTHCO2physeqtethmergedDate)
sample_data(HTHCO2physeqtethmergedDate)$Date
identical(SD,sample_data(HTHCO2physeqtethmergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqtethmergedDate)$Date <- factor(sample_names(HTHCO2physeqtethmergedDate))
#Getting abundance
tethHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqtethmergedDate), "matrix")
View(tethHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
tethHTHCO2_abund_mergedDate_abs_pre<-tethHTHCO2_abund_mergedDate
write.csv(tethHTHCO2_abund_mergedDate, "tethHTHCO2_abund_mergedDate.csv", row.names = F)

###Subset with Verongida
HTHCO2physeqver=subset_taxa(HTHCO2physeqpor,Subclass=="Verongimorpha")
HTHCO2physeqvermergedDate = merge_samples(HTHCO2physeqver, "Date")
SD= merge_samples(sample_data(HTHCO2physeqver), "Date")
print(SD[, "Date"])
print(HTHCO2physeqvermergedDate)
sample_names(HTHCO2physeqvermergedDate)
sample_data(HTHCO2physeqvermergedDate)$Date
identical(SD,sample_data(HTHCO2physeqvermergedDate))
##Merge_samples is a great tool but comes with problem. It changes the values of the metadata and sample names to numerical form. WE have to change the sample_names back to factor by doing the following
sample_data(HTHCO2physeqvermergedDate)$Date <- factor(sample_names(HTHCO2physeqvermergedDate))
#Getting abundance
verHTHCO2_abund_mergedDate=as(otu_table(HTHCO2physeqvermergedDate), "matrix")
View(verHTHCO2_abund_mergedDate)
#incorporating abundance changing it to presence absence data
verHTHCO2_abund_mergedDate_abs_pre<-verHTHCO2_abund_mergedDate
write.csv(verHTHCO2_abund_mergedDate, "verHTHCO2_abund_mergedDate.csv", row.names = F)


#merge all of the csv per sponge group and bring in to R
poriferabyclasorder_HTHCO2_abund_mergedDate<-read.csv("poriferabyclasorder_HTHCO2_abund_mergedDate.csv") #import .csv poriferabyclasorder_Intake_abund_mergedDate.csv
View(poriferabyclasorder_HTHCO2_abund_mergedDate)
poriferabyclasorder_HTHCO2_abund_mergedDate[is.na(poriferabyclasorder_HTHCO2_abund_mergedDate)] <- 0
View(poriferabyclasorder_HTHCO2_abund_mergedDate)
Porif_HTHCO2_class_order<-poriferabyclasorder_HTHCO2_abund_mergedDate
View(Porif_HTHCO2_class_order)
####rarefraction curve###
Porif_HTHCO2_class_order_all<-Porif_HTHCO2_class_order[2:94]
View(Porif_HTHCO2_class_order_all)
curve_Porif_HTHCO2_class_order_all = specaccum(Porif_HTHCO2_class_order_all, method = "rarefaction", 
                                               permutations = 100)
#subset each habitat into its own df
Porif_HTHCO2_class_order%>% filter(Site == "Calcarea") -> Calcarea
names(Calcarea)
Porif_HTHCO2_class_order%>% filter(Site == "Homoscleromorpha") -> Homoscleromorpha
Porif_HTHCO2_class_order%>% filter(Site == "Keratosa") -> Keratosa
Porif_HTHCO2_class_order%>% filter(Site == "Tetractinellida") -> Tetractinellida
Porif_HTHCO2_class_order%>% filter(Site == "Suberitida") -> Suberitida
Porif_HTHCO2_class_order%>% filter(Site == "Poecilosclerida") -> Poecilosclerida
Porif_HTHCO2_class_order%>% filter(Site == "Haplosclerida") -> Haplosclerida
Porif_HTHCO2_class_order%>% filter(Site == "Tethyida") -> Tethyida
Porif_HTHCO2_class_order%>% filter(Site == "Dictyoceratida") -> Dictyoceratida
Porif_HTHCO2_class_order%>% filter(Site == "Dendroceratida") -> Dendroceratida
Porif_HTHCO2_class_order%>% filter(Site == "Verongida") -> Verongida
#species accumulation curve for each habitat using all sponges
curve_HTHCO2_calc = specaccum(Calcarea[, 2:94], method = "rarefaction")
curve_HTHCO2_hom = specaccum(Homoscleromorpha[, 2:94], method = "rarefaction")
curve_HTHCO2_ker = specaccum(Keratosa[, 2:94], method = "rarefaction")
curve_HTHCO2_tetr = specaccum(Tetractinellida[, 2:94], method = "rarefaction")
curve_HTHCO2_sub = specaccum(Suberitida[, 2:94], method = "rarefaction")
curve_HTHCO2_poec = specaccum(Poecilosclerida[, 2:94], method = "rarefaction")
curve_HTHCO2_hap = specaccum(Haplosclerida[, 2:94], method = "rarefaction")
curve_HTHCO2_teth = specaccum(Tethyida[, 2:94], method = "rarefaction")
curve_HTHCO2_ver = specaccum(Verongida[, 2:94], method = "rarefaction")
curve_HTHCO2_Dicty = specaccum(Dictyoceratida[, 2:94], method = "rarefaction")
curve_HTHCO2_Dendr = specaccum(Dendroceratida[, 2:94], method = "rarefaction")
#figure out color palet "display.brewer.pal(n = 8, name = 'RdBu')"
#got the error Error in plot.new() : figure margins too large" so I typed the following command "par(mar=c(5,5,3,1))"" # this fixed the margins and shows axis vallues. 
#plot curve_all first by dates
par(mar=c(5,5,3,1))
plot(curve_HTHCO2_calc, xvar=c("individuals"),ylim=c(0,15),xlim=c(0,800),col="#B2182B", lwd=2, ci.lty=0, ci.col="#F4A582",
     main = "Default: Prettier CI")

#then plot the rest
plot(curve_HTHCO2_hom, add = TRUE,xvar=c("individuals"), col="#2166AC", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI") #col is COLOUR setting, so change it to something else if you 

plot(curve_HTHCO2_ker, add = TRUE,xvar=c("individuals"), col="#E1BE6A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_tetr, add = TRUE,xvar=c("individuals"), col="#40B0A6", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_sub, add = TRUE,xvar=c("individuals"), col="#E66100", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_poec, add = TRUE,xvar=c("individuals"), col="#5D3A9B", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_hap, add = TRUE,xvar=c("individuals"), col="#1AFF1A", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_teth, add = TRUE,xvar=c("individuals"), col="#994F00", lwd=2, ci.lty=0, 
     main = "Default: Prettier CI")



legend(500, 1, legend=c("Calcareous", "Homoscleromorpha"),
       col=c("#B2182B", "#2166AC"), lty=1:1, cex=0.8)

#decided to exclude Dictyo Dendro and Verong since counts were minimal.
```


```{r setup, include=FALSE}
#picking colors 
display.brewer.all(colorblindFriendly=TRUE)

#combining Ambient, HTACO2, ATHCO2 and HTHCO2 abundance in excel... files that have the name porifera"TMThere"_abund_mergedDate.csv
#still need to edit this section

#Bring data of ambient and intake in. 
alltreatmentsrare<-read.csv("porifera_ambient_ATHCO2_HTACO2_HTHCO2_abund_merged.csv") #import .csv
names(alltreatmentsrare)
view(alltreatmentsrare)
alltreatmentsrare_OTUs<-alltreatmentsrare[2:94]
View(alltreatmentsrare_OTUs)
class(alltreatmentsrare_OTUs)
str(alltreatmentsrare_OTUs)
curve_alltreatmentsrare = specaccum(alltreatmentsrare_OTUs, method = "rarefaction", 
                                    permutations = 100)
#subset each habitat into its own df
alltreatmentsrare%>% filter(Site == "ATACO2") -> Ambientabunall

alltreatmentsrare%>% filter(Site == "HTACO2") -> HTACO2abunall

alltreatmentsrare%>% filter(Site == "ATHCO2") -> ATHCO2abunall

alltreatmentsrare%>% filter(Site == "HTHCO2") -> HTHCO2abunall

#species accumulation curve for each habitat using all sponges
curve_HTACO2_abunall = specaccum(HTACO2abunall[, 2:94], method = "rarefaction")
curve_Ambient_abunall = specaccum(Ambientabunall[, 2:94], method = "rarefaction")
curve_ATHCO2_abunall = specaccum(ATHCO2abunall[, 2:94], method = "rarefaction")
curve_HTHCO2_abunall = specaccum(HTHCO2abunall[, 2:94], method = "rarefaction")
curve_HTACO2_abunall


par(mar=c(5,5,3,1))

plot(curve_HTACO2_abunall,ylim=c(30,55), xlim=c(200,2400), xvar=c("individuals"),col="#D55E00", lwd=6, ci.lty=0, ci.col="#E69F00")

#then plot the rest
plot(curve_Ambient_abunall, add = TRUE,xvar=c("individuals"), col="#009E73", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_abunall, add = TRUE,xvar=c("individuals"), col="#56B4E9", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_abunall, add = TRUE,xvar=c("individuals"), col="#999999", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")

##plot curve_all first by dates so 'sites"
par(mar=c(5,5,3,1))

plot(curve_HTACO2_abunall,ylim=c(30,55), xlim=c(1,13), xvar=c("sites"),col="#D55E00", lwd=6, ci.lty=0, ci.col="#E69F00")

#then plot the rest
plot(curve_Ambient_abunall, add = TRUE,xvar=c("sites"), col="#009E73", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_ATHCO2_abunall, add = TRUE,xvar=c("sites"), col="#56B4E9", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")

plot(curve_HTHCO2_abunall, add = TRUE,xvar=c("sites"), col="#999999", lwd=6, ci.lty=0, 
     main = "Default: Prettier CI")


#merging matrices together. Was running into an error with merge but now am using rbind. Issue with the class of data.

#Demospongiae
#rbind to merge dataframes.rbind works with matrices and merge works with dataframes. 
demalltmtmerged<-rbind(demHTHCO2_abund_mergedDate, demATHCO2_abund_mergedDate, demHTACO2_abund_mergedDate, demAmbient_abund_mergedDate)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 
demalltmtmergeddf<-as.data.frame(demalltmtmerged) #Had to change from matrix to dataframe
write_csv(demalltmtmergeddf, "demalltmtmergeddf.csv")
demalltmtmergeddfreimp<-read.csv("demalltmtmergeddf.csv") #import .csv
names(demalltmtmergeddfreimp)
rownames(demalltmtmergeddfreimp)<-NULL #gets rid of rownames
names(demalltmtmergeddfreimp)
demalltmtmergeddfreimp_OTUs<-demalltmtmergeddfreimp[1:66]
names(demalltmtmergeddfreimp_OTUs)
str(demalltmtmergeddfreimp_OTUs)
curve_alltreatmentsraredem = specaccum(demalltmtmergeddfreimp_OTUs, method = "rarefaction", 
                                       permutations = 100) 
#subset each habitat into its own df
demalltmtmergeddfreimp%>% filter(TMT == "ATACO2") -> ATACO2abundemall

demalltmtmergeddfreimp%>% filter(TMT == "HTACO2") -> HTACO2abundemall

demalltmtmergeddfreimp%>% filter(TMT == "ATHCO2") -> ATHCO2abundemall

demalltmtmergeddfreimp%>% filter(TMT == "HTHCO2") -> HTHCO2abundemall

#species accumulation curve for each habitat using all sponges
curve_ATACO2abundem = specaccum(ATACO2abundemall[, 1:66], method = "rarefaction")
curve_HTACO2abundem = specaccum(HTACO2abundemall[, 1:66], method = "rarefaction")
curve_ATHCO2_abundem = specaccum(ATHCO2abundemall[, 1:66], method = "rarefaction")
curve_HTHCO2_abundem = specaccum(HTHCO2abundemall[, 1:66], method = "rarefaction")

#plot curve_all first by dates
par(mar=c(5,5,3,1))
plot(curve_HTACO2abundem, ylim=c(15,40), xlim=c(100,1500), xvar=c("individuals"), col="#D55E00", lwd=4, ci.lty=0, ci.col="#F4A582")


#then plot the rest
plot(curve_ATACO2abundem, add = TRUE,xvar=c("individuals"), col="#009E73", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

plot(curve_ATHCO2_abundem, add = TRUE,xvar=c("individuals"), col="#56B4E9", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

plot(curve_HTHCO2_abundem, add = TRUE,xvar=c("individuals"), col="#999999", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

#Calcarea
#rbind to merge dataframes.rbind works with matrices and merge works with dataframes. 
calcalltmtmerged<-rbind(calcHTHCO2_abund_mergedDate, calcATHCO2_abund_mergedDate, calcHTACO2_abund_mergedDate, calcAmbient_abund_mergedDate)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 

calcalltmtmergeddf<-as.data.frame(calcalltmtmerged) #Had to change from matrix to dataframe
write_csv(calcalltmtmergeddf, "calcalltmtmergeddf.csv")
names(calcalltmtmergeddf)
calcalltmtmergeddfreimp<-read.csv("calcalltmtmergeddf.csv") #import .csv
names(calcalltmtmergeddfreimp)
rownames(calcalltmtmergeddfreimp)<-NULL #gets rid of rownames
names(calcalltmtmergeddfreimp)
calcalltmtmergeddfreimp_OTUs<-calcalltmtmergeddfreimp[1:19]
names(calcalltmtmergeddfreimp_OTUs)
curve_alltreatmentsrarecalc = specaccum(calcalltmtmergeddfreimp_OTUs, method = "rarefaction", 
                                        permutations = 100) 
#subset each habitat into its own df
calcalltmtmergeddfreimp%>% filter(TMT == "ATACO2") -> ATACO2abuncalcall

calcalltmtmergeddfreimp%>% filter(TMT == "HTACO2") -> HTACO2abuncalcall

calcalltmtmergeddfreimp%>% filter(TMT == "ATHCO2") -> ATHCO2abuncalcall

calcalltmtmergeddfreimp%>% filter(TMT == "HTHCO2") -> HTHCO2abuncalcall

#species accumulation curve for each habitat using all sponges
curve_ATACO2abuncalc = specaccum(ATACO2abuncalcall[, 1:19], method = "rarefaction")
curve_HTACO2abuncalc = specaccum(HTACO2abuncalcall[, 1:19], method = "rarefaction")
curve_ATHCO2_abuncalc = specaccum(ATHCO2abuncalcall[, 1:19], method = "rarefaction")
curve_HTHCO2_abuncalc = specaccum(HTHCO2abuncalcall[, 1:19], method = "rarefaction")

#plot curve_all first by dates
par(mar=c(5,5,3,1))

plot(curve_ATACO2abuncalc,ylim=c(10,15), xlim=c(50,800), xvar=c("individuals"),col="#009E73", lwd=4, ci.lty=0, ci.col="#F4A582")

#then plot the rest
plot(curve_ATHCO2_abuncalc, add = TRUE,xvar=c("individuals"), col="#56B4E9", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

plot(curve_HTACO2abuncalc, add = TRUE,xvar=c("individuals"), col="#D55E00", lwd=4, ci.lty=0) 


plot(curve_HTHCO2_abuncalc, add = TRUE,xvar=c("individuals"), col="#999999", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

#Homoscleromorpha
#Homoscleromorpha
#rbind to merge dataframes.rbind works with matrices and merge works with dataframes. 
homoalltmtmerged<-rbind(homoHTHCO2_abund_mergedDate, homoATHCO2_abund_mergedDate, homoHTACO2_abund_mergedDate, homoAmbient_abund_mergedDate)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 

homoalltmtmergeddf<-as.data.frame(homoalltmtmerged) #Had to change from matrix to dataframe
write_csv(homoalltmtmergeddf, "homoalltmtmergeddf.csv")
names(calcalltmtmergeddf)
homoalltmtmergeddfreimp<-read.csv("homoalltmtmergeddf.csv") #import .csv
names(homoalltmtmergeddfreimp)
rownames(homoalltmtmergeddfreimp)<-NULL #gets rid of rownames
names(homoalltmtmergeddfreimp)
homoalltmtmergeddfreimp_OTUs<-homoalltmtmergeddfreimp[1:8]
names(homoalltmtmergeddfreimp_OTUs)
curve_alltreatmentsrarehomo = specaccum(homoalltmtmergeddfreimp_OTUs, method = "rarefaction", 
                                        permutations = 100) 
#subset each habitat into its own df
homoalltmtmergeddfreimp%>% filter(TMT == "ATACO2") -> ATACO2abunhomoall

homoalltmtmergeddfreimp%>% filter(TMT == "HTACO2") -> HTACO2abunhomoall

homoalltmtmergeddfreimp%>% filter(TMT == "ATHCO2") -> ATHCO2abunhomoall

homoalltmtmergeddfreimp%>% filter(TMT == "HTHCO2") -> HTHCO2abunhomoall

#species accumulation curve for each habitat using all sponges
curve_ATACO2abunhomo = specaccum(ATACO2abunhomoall[, 1:8], method = "rarefaction")
curve_HTACO2abunhomo = specaccum(HTACO2abunhomoall[, 1:8], method = "rarefaction")
curve_ATHCO2_abunhomo = specaccum(ATHCO2abunhomoall[, 1:8], method = "rarefaction")
curve_HTHCO2_abunhomo = specaccum(HTHCO2abunhomoall[, 1:8], method = "rarefaction")

#plot curve_all first by dates
par(mar=c(5,5,3,1))

plot(curve_ATACO2abunhomo,ylim=c(0,7), xlim=c(00,200), xvar=c("individuals"),col="#009E73", lwd=4, ci.lty=0, ci.col="#F4A582")

#then plot the rest
plot(curve_HTACO2abunhomo, add = TRUE,xvar=c("individuals"), col="#D55E00", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

plot(curve_ATHCO2_abunhomo, add = TRUE,xvar=c("individuals"), col="#56B4E9", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 

plot(curve_HTHCO2_abunhomo, add = TRUE,xvar=c("individuals"), col="#999999", lwd=4, ci.lty=0) #col is COLOUR setting, so change it to something else if you 




#Richness box blots for all sponges
poralltmtmerged_richness<-rbind(Richness_HTHCO2physeqpormergedunit, Richness_ATHCO2physeqpormergedunit, Richness_HTACO2physeqpormergedunit, Richness_Ambientphyseqpormergedunit)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 
poralltmtmerged_richness_df<-as.data.frame(poralltmtmerged_richness) #Had to change from matrix to dataframe
write_csv(poralltmtmerged_richness_df, "poralltmtmerged_richness_df.csv")
poralltmtmergeddfrichness<-read.csv("poralltmtmerged_richness_df.csv") #import .csv
poralltmtmergeddfrichness
poralltmtmergeddfrichness<-cbind(poralltmtmergeddfrichness, Group='Porifera')
porcumrichnessplot<-ggplot(poralltmtmergeddfrichness, aes(x = factor(Observed, level = c('ATACO2', 'HTACO2', 'ATHCO2','HTHCO2')), y=Observed,
                                                          fill=TMT)) + theme_classic()+
  geom_boxplot(outlier.colour = "red")+stat_summary(fun=mean, geom="point", position=position_dodge(width=0.75), 
                                                    color="black", size=2)+
  labs(y = "Observed richness", x = element_blank())+scale_fill_manual("legend", values=c("ATACO2" = "#009E73","HTACO2" ="#D55E00","ATHCO2" = "#56B4E9", "HTHCO2" = "#999999"))