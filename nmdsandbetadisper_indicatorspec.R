#preparing data for NMDS analysis
#separating by treatment and then merging to create an OTU table that is usable
#ambient
Ambientphyseqpormergedunit_nmds = merge_samples(Ambientphyseqpor, "UNIT")
SD= merge_samples(sample_data(Ambientphyseqpor), "UNIT")
print(SD[, "UNIT"])
print(Ambientphyseqpormergedunit_nmds)
sample_names(Ambientphyseqpormergedunit_nmds)
sample_data(Ambientphyseqpormergedunit_nmds)$UNIT
identical(SD,sample_data(Ambientphyseqpormergedunit_nmds))
Ambientphyseqpormergedunit_nmds_otu=as(otu_table(Ambientphyseqpormergedunit_nmds), "matrix")
Ambientphyseqpormergedunit_nmds_otu
Ambientphyseqpormergedunit_nmds_otu<-cbind(Ambientphyseqpormergedunit_nmds_otu, TMT='ATACO2')
row.names(Ambientphyseqpormergedunit_nmds_otu)=sample_names(Ambientphyseqpormergedunit_nmds)
Ambientphyseqpormergedunit_nmds_otu<-dfRowName(Ambientphyseqpormergedunit_nmds_otu, name ="UNIT")
Ambientphyseqpormergedunit_nmds_otu

#HTACO2
HTACO2physeqpormergedunit_nmds = merge_samples(HTACO2physeqpor, "UNIT")
SD= merge_samples(sample_data(HTACO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(HTACO2physeqpormergedunit_nmds)
sample_names(HTACO2physeqpormergedunit_nmds)
sample_data(HTACO2physeqpormergedunit_nmds)$UNIT
identical(SD,sample_data(HTACO2physeqpormergedunit_nmds))
HTACO2physeqpormergedunit_nmds_otu=as(otu_table(HTACO2physeqpormergedunit_nmds), "matrix")
HTACO2physeqpormergedunit_nmds_otu
HTACO2physeqpormergedunit_nmds_otu<-cbind(HTACO2physeqpormergedunit_nmds_otu, TMT='HTACO2')
row.names(HTACO2physeqpormergedunit_nmds_otu)=sample_names(HTACO2physeqpormergedunit_nmds)
HTACO2physeqpormergedunit_nmds_otu<-dfRowName(HTACO2physeqpormergedunit_nmds_otu, name ="UNIT")
HTACO2physeqpormergedunit_nmds_otu

#ATHCO2
ATHCO2physeqpormergedunit_nmds = merge_samples(ATHCO2physeqpor, "UNIT")
SD= merge_samples(sample_data(ATHCO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(ATHCO2physeqpormergedunit_nmds)
sample_names(ATHCO2physeqpormergedunit_nmds)
sample_data(ATHCO2physeqpormergedunit_nmds)$UNIT
identical(SD,sample_data(ATHCO2physeqpormergedunit_nmds))
ATHCO2physeqpormergedunit_nmds_otu=as(otu_table(ATHCO2physeqpormergedunit_nmds), "matrix")
ATHCO2physeqpormergedunit_nmds_otu
ATHCO2physeqpormergedunit_nmds_otu<-cbind(ATHCO2physeqpormergedunit_nmds_otu, TMT='ATHCO2')
row.names(ATHCO2physeqpormergedunit_nmds_otu)=sample_names(ATHCO2physeqpormergedunit_nmds)
ATHCO2physeqpormergedunit_nmds_otu<-dfRowName(ATHCO2physeqpormergedunit_nmds_otu, name ="UNIT")
ATHCO2physeqpormergedunit_nmds_otu


#HTHCO2
HTHCO2physeqpormergedunit_nmds = merge_samples(HTHCO2physeqpor, "UNIT")
SD= merge_samples(sample_data(HTHCO2physeqpor), "UNIT")
print(SD[, "UNIT"])
print(HTHCO2physeqpormergedunit_nmds)
sample_names(HTHCO2physeqpormergedunit_nmds)
sample_data(HTHCO2physeqpormergedunit_nmds)$UNIT
identical(SD,sample_data(HTHCO2physeqpormergedunit_nmds))
HTHCO2physeqpormergedunit_nmds_otu=as(otu_table(HTHCO2physeqpormergedunit_nmds), "matrix")
HTHCO2physeqpormergedunit_nmds_otu
HTHCO2physeqpormergedunit_nmds_otu<-cbind(HTHCO2physeqpormergedunit_nmds_otu, TMT='HTHCO2')
row.names(HTHCO2physeqpormergedunit_nmds_otu)=sample_names(HTHCO2physeqpormergedunit_nmds)
HTHCO2physeqpormergedunit_nmds_otu<-dfRowName(HTHCO2physeqpormergedunit_nmds_otu, name ="UNIT")

#merged all tmts
alltmtmergedunitOTUmds<-rbind(HTHCO2physeqpormergedunit_nmds_otu, Ambientphyseqpormergedunit_nmds_otu, ATHCO2physeqpormergedunit_nmds_otu, HTACO2physeqpormergedunit_nmds_otu)
#as dataframe
alltmtmergedunitOTUmdsdf<-as.data.frame(alltmtmergedunitOTUmds)
view(alltmtmergedunitOTUmds)
alltmtmergedunitOTUmdsdf= merge(HTcol_df, alltmtmergedunitOTUmds, by="UNIT")#merging HT
view(alltmtmergedunitOTUmdsdf)
write_csv(alltmtmergedunitOTUmdsdf, "alltmtmergedunitOTUmdsdf.csv")
alltmtmergedunitOTUmdsdf<-read.csv("alltmtmergedunitOTUmdsdf.csv") #import .csv
names(alltmtmergedunitOTUmdsdf)
alltmtmergedunitOTUmdsdf$TMT.y<-NULL
names(alltmtmergedunitOTUmdsdf)
colnames(alltmtmergedunitOTUmdsdf)[3] ="TMT" #change column name back to TMT
names(alltmtmergedunitOTUmdsdf)
metamergedunit<-alltmtmergedunitOTUmdsdf[,-(4:96)] # made metadata file by deleting OTUs from alltmtmergedOTUmdsdf
names(metamergedunit)
view(metamergedunit)
str(metamergedunit)



#For vegan porifera
view(alltmtmergedunitOTUmdsdf)
metavallTMT=metamergedunit
OTUallTMT=alltmtmergedunitOTUmdsdf #
Tax=OTUClass_df.mat
View(Tax)
names(OTUallTMT)
OTUallTMT<-OTUallTMT[,-1]
names(OTUallTMT)
OTUallTMT<-OTUallTMT[,-(2)]
names(OTUallTMT)
OTUallTMT<-OTUallTMT[,-(1)]
names(OTUallTMT)
View(OTUallTMT)
#Remove OTUs that dont exits in OTUvallTMT
OTUallTMTc<-OTUallTMT[,-(which(colSums(OTUallTMT)==0))]
View(OTUallTMTc)
#remove all the OTUs that don't occur in OTU data set. 
TaxallTMT.clean = Tax[row.names(Tax) %in% colnames(OTUallTMT),]
View(TaxallTMT.clean)
#Order the data
OTUallTMTco = OTUallTMTc[order(row.names(OTUallTMT)),]
View(OTUallTMTco)
metaallTMTo=metavallTMT[order(row.names(metavallTMT)),]
View(metaallTMTo)
str(metaallTMTo)
#mutating to include pco2
metaallTMTo <- metaallTMTo %>% 
  mutate(TMT2 = TMT)
names(metaallTMTo)
colnames(metaallTMTo)[4] <- "pco2" #change colnames of   TMT2 to temp
names(metaallTMTo)
metaallTMTo
metaallTMTo$pco2 <- revalue(metaallTMTo$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
names(metaallTMTo)

#mutating to iclude temp
metaallTMTo <- metaallTMTo %>% 
  mutate(TMT2 = TMT)
names(metaallTMTo)
colnames(metaallTMTo)[5] <- "temp" #change colnames of   TMT2 to temp
names(metaallTMTo)
view(metaallTMTo)
class(metaallTMTo)
metaallTMTo$temp <- revalue(metaallTMTo$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
view(metaallTMTo)
names(metaallTMTo)
str(metaallTMTo)
#for relative abundance stack plot
allTMTporrelabu<-cbind(OTUallTMTco, metaallTMTo)
head(allTMTporrelabu)
#plot relative abundance
relabunplot<-ggplot(allTMTporrelabu, aes(x = TMT, y = Abundance, fill = TaxallTMT.clean$Class)) + 
  geom_bar(stat = "identity")
relabunplot

#Transform the data sqrt
view(OTUallTMTco)
OTUallTMTcosq<-sqrt(OTUallTMTco)
View(OTUallTMTcosq)
names(OTUallTMTcosq)
#for pairwise analysis cbindg OTUallTMTcosq and 
allTMTpairwise<-cbind(OTUallTMTcosq, metaallTMTo)
names(allTMTpairwise)
#Bray curtis distances
BC.nmds_allTMT<-metaMDS(OTUallTMTcosq, distance="bray",autotransform = FALSE) 
BC.nmds_allTMT#Checking stress value: Stress:       0.2313  in this case
# Shepards test/goodness of fit
goodness(BC.nmds_allTMT) # Produces a results of test statistics for goodness of fit for each point

stressplot(BC.nmds_allTMT) # Produces a Shepards diagram
plot(BC.nmds_allTMT) # displays habitats and species
view(metaallTMTo)
metaallTMTo$Header<-as.factor(metaallTMTo$Header)

str(metaallTMTo)
#running permanova with nesting as in https://stackoverflow.com/questions/77884788/nested-permanova-in-r
permanovaTMTwithHTnested1 <-adonis2(formula = vegdist(OTUallTMTcosq, distance = "bray") ~ TMT, data = metaallTMTo, strata = metaallTMTo$Header, permutations = 1000)
permanovaTMTwithHTnested1
permanovaTMTwithHTnested2 <-adonis2(formula = vegdist(OTUallTMTcosq, distance = "bray") ~ TMT*Header, data = metaallTMTo, , permutations = 1000, by = 'margin')
permanovaTMTwithHTnested2
permanovaTMTwithHTnested3 <-adonis2(formula = vegdist(OTUallTMTcosq, distance = "bray") ~ TMT/Header, data = metaallTMTo, , permutations = 1000, by = 'terms')
permanovaTMTwithHTnested3

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 1000

adonis2(formula = vegdist(OTUallTMTcosq, distance = "bray") ~ TMT/Header, data = metaallTMTo, permutations = 1000, by = "terms")
Df SumOfSqs      R2      F   Pr(>F)    
TMT         3  0.60906 0.23195 2.0994 0.000999 ***
  TMT:Header  4  0.46949 0.17880 1.2137 0.080919 .  
Residual   16  1.54727 0.58925                    
Total      23  2.62583 1.00000                    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

print(permanovaTMTwithHTnested3)

#This is it! export to excel! using dplyr
permanovaTMTwithHTnested3_df <- tidy(permanovaTMTwithHTnested3) %>% 
  mutate(data.name = permanovaTMTwithHTnested3$data.name, 
         dataset = "metaallTMTo")
permanovaTMTwithHTnested3_df
write_csv(permanovaTMTwithHTnested3_df, "permanovaTMTwithHTnested3_df.csv")

names(permanovaTMTwithHTnested3) <- colnames(metaallTMTo)
permanovaTMTwithHTnested3
permanovaTMTwithHTnested3 <- do.call(rbind, permanovaTMTwithHTnested3)
permanovaTMTwithHTnested3

# Pairwise comparisons between Management levels
metaallTMTo$TMT=as.factor(metaallTMTo$TMT)
# Pairwise comparisons between Management levels
pwnmdsTMT <- adonis_pairwise(x = metaallTMTo, dd = vegdist(OTUallTMTcosq), group.var = "TMT")
pwnmdsTMT

$Adonis.tab
Comparison        R2        F   df     p p.adj
1 ATACO2.ATHCO2 0.1344676 1.553583 1;10 0.091 0.091
2 ATACO2.HTACO2 0.1670089 2.004930 1;10 0.021 0.030
3 ATACO2.HTHCO2 0.2034400 2.553983 1;10 0.005 0.014
4 ATHCO2.HTACO2 0.1405607 1.635493 1;10 0.025 0.030
5 ATHCO2.HTHCO2 0.1988709 2.482383 1;10 0.001 0.006
6 HTACO2.HTHCO2 0.1604488 1.911125 1;10 0.007 0.014

$Betadisper.tab
Comparison          F   df     p  p.adj
1 ATACO2.ATHCO2 2.05458030 1;10 0.177 0.2655
2 ATACO2.HTACO2 0.27525940 1;10 0.619 0.7428
3 ATACO2.HTHCO2 2.34530147 1;10 0.155 0.2655
4 ATHCO2.HTACO2 2.75119228 1;10 0.117 0.2655
5 ATHCO2.HTHCO2 0.01176429 1;10 0.870 0.8700
6 HTACO2.HTHCO2 2.96608115 1;10 0.090 0.2655


#following https://www.rpubs.com/RGrieger/545184
plot(BC.nmds_allTMT, type = "n") #displays empty ordination space
points(BC.nmds_allTMT, display = "sites", pch = c(16, 8, 4, 15) [as.numeric(metaallTMTo$TMT)], col = c("blue", "orange", "green", "black") [as.numeric(metaallTMTo$TMT)])# displays site points where symbols (pch) are different management options and colour (col) are different land uses
legend("top", legend = c(levels(metaallTMTo$TMT)), pch = c(16, 8, 4, 15), col = c("blue", "orange", "green", "black"), bty = "n", cex = 1)
legend("bottom", "stress = 0.23", bty = "n", cex = 1) # displays legend text of stress value 


#Jaccard curtis distances
JC.nmds_allTMT<-metaMDS(OTUallTMTcosq, distance="jaccard",autotransform = FALSE) 
JC.nmds_allTMT#Checking stress value: Stress:     0.23in this case
plot(JC.nmds_allTMT) # displays habitats and species

#following https://www.rpubs.com/RGrieger/545184
plot(JC.nmds_allTMT, type = "n") #displays empty ordination space
points(JC.nmds_allTMT, display = "sites", pch = c(16, 8, 4, 15) [as.numeric(metaallTMTo$TMT)], col = c("blue", "orange", "green", "black") [as.numeric(metaallTMTo$TMT)]) # displays site points where symbols (pch) are different management options and colour (col) are different land uses
legend("top", legend = c(levels(metaallTMTo$TMT)), pch = c(16, 8, 4, 15), col = c("blue", "orange", "green", "black"), bty = "n", cex = 1) # displays symbol and colour legend
legend("bottom", "stress = 0.23", bty = "n", cex = 1) # displays legend text of stress value 
#investigate the species driving distribution
allTMT.spp.fit <- envfit(JC.nmds_allTMT, OTUallTMTcosq, permutations = 999)
head(allTMT.spp.fit)

ordiplot(JC.nmds_allTMT, type = "n", main = "intrinsic species")

orditorp(JC.nmds_allTMT, display = "sites", labels = F, pch = c(16, 8) [as.numeric(metaallTMTo$DATE)], cex = 1)

plot(allTMT.spp.fit, p.max = 0.01, col = "black", cex = 0.7) # change the significance level of species shown with p.max

#PERMANOVA on temp

#PERMANOVA for braycurtis distances but first need to add "HT" column to the "metallTMTo" metadata file. 
#importing
HTcol<-read.csv("HT.csv", header = T, strip.white = T, na.strings = "")
#making HTcol a dataframe
HTcol_df<-as.data.frame(HTcol)
view(HTcol_df)
#merging df based on UNIT
metaallTMTo_HT= merge(HTcol_df, metaallTMTo, by="UNIT")
view(metaallTMTo_HT)
metaallTMTo_HT$Header<-as.factor(metaallTMTo_HT$Header)
metaallTMTo_HT$TMT<-as.factor(metaallTMTo_HT$TMT)
metaallTMTo_HT$pco2<-as.factor(metaallTMTo_HT$pco2)
metaallTMTo_HT$temp<-as.factor(metaallTMTo_HT$temp)
str(metaallTMTo_HT) #make sure that Header is a factor!


#with nesting
BC.nmds_allTMT_dist_adonis_temp_pco2_HT<-adonis2(formula = BC.nmds_allTMT_dist ~ temp*pco2/factor(Header), data = metaallTMTo_HT, permutations = 1000)
BC.nmds_allTMT_dist_adonis_temp_pco2_HT

#PERMANOVA on temp and pco2 with nesting
permanovapco2andTwithHTnested3 <-adonis2(formula = vegdist(OTUallTMTcosq, distance = "bray") ~ temp*pco2/Header, data = metaallTMTo, permutations = 1000, by = 'terms')
permanovapco2andTwithHTnested3

#This is it! export to excel! using dplyr
permanovapco2andTwithHTnested3_df <- tidy(permanovapco2andTwithHTnested3) %>% 
  mutate(data.name = permanovapco2andTwithHTnested3$data.name, 
         dataset = "metaallTMTo")
permanovapco2andTwithHTnested3_df
write_csv(permanovapco2andTwithHTnested3_df, "permanovapco2andTwithHTnested3_df.csv")

# A tibble: 6 × 7
term                df SumOfSqs     R2 statistic   p.value dataset    
<chr>            <dbl>    <dbl>  <dbl>     <dbl>     <dbl> <chr>      
  1 temp                 1    0.259 0.0986      2.68  0.000999 metaallTMTo
2 pco2                 1    0.162 0.0616      1.67  0.0180   metaallTMTo
3 temp:pco2            1    0.188 0.0718      1.95  0.00300  metaallTMTo
4 temp:pco2:Header     4    0.469 0.179       1.21  0.111    metaallTMTo
5 Residual            16    1.55  0.589      NA    NA        metaallTMTo
6 Total               23    2.63  1          NA    NA        metaallTMTo


#betadisper for temperature
# Pairwise comparisons between Management levels
metaallTMTo$temp=as.factor(metaallTMTo$temp)
# Pairwise comparisons between Management levels
pwnmdstemp <- adonis_pairwise(x = metaallTMTo, dd = vegdist(OTUallTMTcosq, distance = "bray"), group.var = "temp")
pwnmdstemp
$Adonis.tab
Comparison         R2        F   df     p p.adj
1 Future.Present 0.09855813 2.405345 1;22 0.002 0.002

$Betadisper.tab
Comparison         F   df     p p.adj
1 Future.Present 0.4263081 1;22 0.539 0.539

#betadisper for pCO2
metaallTMTo$pco2=as.factor(metaallTMTo$pco2)
# Pairwise comparisons between Management levels
pwnmdspco2 <- adonis_pairwise(x = metaallTMTo, dd = vegdist(OTUallTMTcosq, distance = "bray"), group.var = "pco2")
pwnmdspco2
$Adonis.tab
Comparison         R2        F   df     p p.adj
1 Future.Present 0.06163412 1.445013 1;22 0.074 0.074

$Betadisper.tab
Comparison        F   df     p p.adj
1 Future.Present 3.028064 1;22 0.095 0.095


#PERMANOVA for Jaccard distances 
JC.nmds_allTMT_dist=vegdist(OTUallTMTcosq, distance="jaccard")
JC.nmds_allTMT_distadonis<-adonis2(formula = JC.nmds_allTMT_dist ~ TMT, data = metaallTMTo, permutations = 1000)
JC.nmds_allTMT_distadonis

Df SumOfSqs      R2      F Pr(>F)
TMT       3   0.6512 0.22267 1.9097 0.000999 ***
  Residual 20   2.2733 0.77733                    
Total    23   2.9245 1.00000                    
--- 
  pairwise.adonis2(JC.nmds_allTMT_dist ~ TMT, data = metaallTMTo)

adonis_pairwise(JC.nmds_allTMT_dist ~ TMT, data = metaallTMTo)
# Compare all Management levels
adonis2(JC.nmds_allTMT_dist ~ TMT, data = metaallTMTo)

# Pairwise comparisons between Management levels
tst <- adonis_pairwise(x = dune.env, dd = vegdist(dune), group.var = "Management")
tst$Adonis.tab
tst$Betadisper.tab

# Pairwise comparisons between Management levels
tst <- adonis_pairwise(x = metaallTMTo, dd = vegdist(JC.nmds_allTMT_dist), group.var = "TMT")
tst$Adonis.tab
tst$Betadisper.tab






#To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. You can do this by calling the scores of you mds.
allTMT_nmds_sites.scrs <- as.data.frame(scores(allTMT_nmds, display = "sites")) #save NMDS results into dataframe
allTMT_site.scrs <- cbind(allTMT_nmds_sites.scrs, TMT = metaallTMTo$TMT) #add grouping variable "Management" to dataframe
head(allTMT_site.scrs)
allTMT_site.scrs
#copying the same column twice to make a factorial design
allTMT_site.scrs <- allTMT_site.scrs %>% mutate (TMT2 = TMT) 
names(allTMT_site.scrs)
#change colnames of   TMT2 to temp
colnames(allTMT_site.scrs)[4] <- "temp" #change colnames of   TMT2 to temp
allTMT_site.scrs$temp <- revalue(allTMT_site.scrs$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(allTMT_site.scrs)
allTMT_site.scrs <- allTMT_site.scrs %>% 
  mutate(TMT2 = TMT)
names(allTMT_site.scrs)
colnames(allTMT_site.scrs)[5] <- "pco2" #change colnames of   TMT2 to temp
names(allTMT_site.scrs)
allTMT_site.scrs
allTMT_site.scrs$pco2 <- revalue(allTMT_site.scrs$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
allTMT_site.scrs
View(allTMT_site.scrs)

#plot in ggplot
allTMT_nmds <- metaMDS(OTUallTMTcosq, distance = "bray", autotransform = F)
allTMT.envfit <- envfit(allTMT_nmds, metaallTMTo, permutations = 999)
allTMT.spp.fit<- envfit(allTMT_nmds, OTUallTMTcosq, permutations = 999) # this fits species vectors

#To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. You can do this by calling the scores of you mds.
allTMT_nmds_sites.scrs <- as.data.frame(scores(allTMT_nmds, display = "sites")) #save NMDS results into dataframe
allTMT_site.scrs <- cbind(allTMT_nmds_sites.scrs, TMT = metaallTMTo$TMT) #add grouping variable "Management" to dataframe
head(allTMT_site.scrs)
allTMT_site.scrs

#A new dataset containing species data also needs to be made to look at species vectors.This is not necessary if you don’t want to show the species on the final graph.
allTMT.spp.scrs<- as.data.frame(scores(allTMT.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
allTMT.spp.scrs <- cbind(allTMT.spp.scrs, Species = rownames(allTMT.spp.scrs)) #add species names to dataframe
allTMT.spp.scrs <- cbind(allTMT.spp.scrs, pval = allTMT.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant

#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
allTMT_sig.spp.scrs <- subset(allTMT.spp.scrs, pval<=0.01) #subset data to show species significant at 0.01
head(allTMT_sig.spp.scrs)

#To show environmental extrinsic variables another datasheet needs to be created
metaallTMTo.scores.allTMT <- as.data.frame(scores(allTMT.envfit, display = "vectors")) #extracts relevant scores from envifit
metaallTMTo.scores.allTMT <- cbind(metaallTMTo.scores.allTMT, metavo.variables = rownames(metaallTMTo.scores.allTMT)) #and then gives them their names
metaallTMTo.scores.allTMT <- cbind(metaallTMTo.scores.allTMT , pval = allTMT.envfit$vectors$pvals) # add pvalues to dataframe
allTMT_sig.env.scrs <- subset(metaallTMTo.scores.allTMT, pval<=0.01) #subset data to show variables significant at 0.05

head(metaallTMTo.scores.allTMT)
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!

nmds.plot.metaallTMTo <- ggplot(allTMT_site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(TMT), size = 3))+geom_convexhull(alpha = 0.5, aes(fill=TMT))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic() + scale_color_manual(values=c("#56B4E9", "#FFE93F", "#D55E00", "#999999"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "TMT")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
nmds.plot.metaallTMTo                         
#Significant Species
nmds.plot.metaallTMTowSP<-nmds.plot.metaallTMTo +geom_segment(data = allTMT_sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = allTMT_sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label=Species), cex = 3, direction = "both", segment.size = 0.25)+labs(title = "Ordination with species vectors") #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
nmds.plot.metaallTMTowSP

#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!
#Treatment as factor
nmds.plot.allTMT <- ggplot(allTMT_site.scrs, aes(x=NMDS1, y=NMDS2, col =TMT))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(TMT), size = 3))+geom_convexhull(alpha = 0, aes(fill=TMT))+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic() +
  scale_color_manual(values=c("ATACO2" = "#56B4E9","ATHCO2" = "#FFE93F","HTACO2" ="#D55E00", "HTHCO2" = "#999999"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
nmds.plot.allTMT 

nmds.plot.metaallTMTowSP<-nmds.plot.allTMT+geom_segment(data = allTMT_sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = allTMT_sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label=Species), cex = 3, direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
labs(title = "Ordination with species vectors")

nmds.plot.metaallTMTowSP

#pCO2
allTMT_site.scrs
nmds.plot.allpco2 <- ggplot(allTMT_site.scrs, aes(x=NMDS1, y=NMDS2, col =pco2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(pco2), size = 2))+geom_convexhull(alpha = 0, aes(fill=pco2))+
  coord_fixed()+
  theme_classic() +
  scale_color_manual(values=c("Future" = "grey67","Present" ="navajowhite2"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
nmds.plot.allpco2

#temp
nmds.plot.alltemp <- ggplot(allTMT_site.scrs, aes(x=NMDS1, y=NMDS2, col =temp))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(temp), size = 2))+geom_convexhull(alpha = 0, aes(fill=temp))+
  coord_fixed()+
  theme_classic() +
  scale_color_manual(values=c("Future" = "grey67","Present" ="navajowhite2"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
nmds.plot.alltemp


p<-ggplot(allTMT_site.scrs,aes(x = NMDS1, y = NMDS2,col = TMT)) + scale_color_manual(values=c("ATACO2" = "#56B4E9","HTACO2" ="#D55E00","ATHCO2" = "#FFE93F", "HTHCO2" = "#999999"))+
  geom_point() +
  theme_minimal()+geom_polygon(alpha = 0.3, aes(fill=TMT))

pp2<-p+geom_convexhull(alpha = 0.3,aes(fill = TMT))

pp2
p
p2

#plot in ggplot with indicator species
#treatments
allnmds_indsp <- metaMDS(OTUallTMTcosq, distance = "bray", autotransform = F)
allnmds_indsp.envfit <- envfit(allnmds_indsp , metaallTMTo, permutations = 999) # this fits environmental vectors
allnmds_indsp.spp.fit <- envfit(allnmds_indsp, OTUallTMTcosq, permutations = 999) # this fits species vectors

#To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. You can do this by calling the scores of you mds.
allnmds_indsp_site.scrs <- as.data.frame(scores(allnmds_indsp, display = "sites")) #save NMDS results into dataframe
allnmds_indsp_site.scrs <- cbind(allnmds_indsp_site.scrs, TMT = metaallTMTo$TMT) #add grouping variable "Management" to dataframe
head(allnmds_indsp_site.scrs)

#A new dataset containing species data also needs to be made to look at species vectors.This is not necessary if you don’t want to show the species on the final graph.
allnmds_indsp_spp.scrs <- as.data.frame(scores(allnmds_indsp.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
allnmds_indsp_spp.scrs <- cbind(allnmds_indsp_spp.scrs, Species = rownames(allnmds_indsp_spp.scrs)) #add species names to dataframe
allnmds_indsp_spp.scrs <- cbind(allnmds_indsp_spp.scrs, pval = allnmds_indsp.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant

#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
allnmds_indsp_sig.spp.scrs <- subset(allnmds_indsp_spp.scrs , pval<=0.01) #subset data to show species significant at 0.01
head(allnmds_indsp_sig.spp.scrs)


#To show environmental extrinsic variables another datasheet needs to be created
metavo.scores.allnmds_indsp <- as.data.frame(scores(allnmds_indsp.envfit, display = "vectors")) #extracts relevant scores from envifit
metavo.scores.allnmds_indsp <- cbind(metavo.scores.allnmds_indsp, metavo.variables = rownames(metavo.scores.allnmds_indsp)) #and then gives them their names

metavo.scores.allnmds_indsp <- cbind(metavo.scores.allnmds_indsp , pval = allnmds_indsp.envfit$vectors$pvals) # add pvalues to dataframe
allnmds_indsp_sig.env.scrs <- subset(metavo.scores.allnmds_indsp, pval<=0.01) #subset data to show variables significant at 0.05
head(metavo.scores.allnmds_indsp)
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!


nmds.plot.allnmds_indsp  <- ggplot(allnmds_indsp_site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(TMT), size = 2))+geom_convexhull(alpha = 0.3, aes(fill=TMT))+#adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic() + 
  scale_color_manual(values=c("ATACO2" = "#009E73","HTACO2" ="#D55E00","ATHCO2" = "#56B4E9", "HTHCO2" = "#999999"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot.allnmds_indsp 

nmds.plot.allnmds_indsp_vectors<-nmds.plot.allnmds_indsp+ geom_segment(data = allnmds_indsp_sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = allnmds_indsp_sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label=Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")
nmds.plot.allnmds_indsp_vectors

#By temperature
allnmds_indsp <- metaMDS(OTUallTMTcosq, distance = "bray", autotransform = F)
allnmds_indsp.envfit <- envfit(allnmds_indsp , metaallTMTo, permutations = 999) # this fits environmental vectors
allnmds_indsp.spp.fit <- envfit(allnmds_indsp, OTUallTMTcosq, permutations = 999) # this fits species vectors

#To plot the output from the mds using ggplot a new datasheet needs to be created which contains the x,y points for each site. You can do this by calling the scores of you mds.
allnmds_indsp_site.scrs_temp <- as.data.frame(scores(allnmds_indsp, display = "sites")) #save NMDS results into dataframe
allnmds_indsp_site.scrs_temp <- cbind(allnmds_indsp_site.scrs, temp = metaallTMTo$temp) #add grouping variable "Management" to dataframe
head(allnmds_indsp_site.scrs_temp)

#A new dataset containing species data also needs to be made to look at species vectors.This is not necessary if you don’t want to show the species on the final graph.
allnmds_indsp_spp.scrs <- as.data.frame(scores(allnmds_indsp.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
allnmds_indsp_spp.scrs <- cbind(allnmds_indsp_spp.scrs, Species = rownames(allnmds_indsp_spp.scrs)) #add species names to dataframe
allnmds_indsp_spp.scrs <- cbind(allnmds_indsp_spp.scrs, pval = allnmds_indsp.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant

#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
allnmds_indsp_sig.spp.scrs <- subset(allnmds_indsp_spp.scrs , pval<=0.01) #subset data to show species significant at 0.05
head(allnmds_indsp_sig.spp.scrs)


#To show environmental extrinsic variables another datasheet needs to be created
metavo.scores.allnmds_indsp <- as.data.frame(scores(allnmds_indsp.envfit, display = "vectors")) #extracts relevant scores from envifit
metavo.scores.allnmds_indsp <- cbind(metavo.scores.allnmds_indsp, metavo.variables = rownames(metavo.scores.allnmds_indsp)) #and then gives them their names

metavo.scores.allnmds_indsp <- cbind(metavo.scores.allnmds_indsp , pval = allnmds_indsp.envfit$vectors$pvals) # add pvalues to dataframe
allnmds_indsp_sig.env.scrs <- subset(metavo.scores.allnmds_indsp, pval<=0.01) #subset data to show variables significant at 0.05
head(metavo.scores.allnmds_indsp)
#Now we have the relevant information for plotting the ordination in ggplot! Lets get plotting!


nmds.plot.allnmds_indsp_temp  <- ggplot(allnmds_indsp_site.scrs_temp, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(temp), size = 2))+geom_convexhull(alpha = 0.3, aes(fill=temp))+#adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic() + 
  scale_color_manual(values=c("Future" = "grey67","Present" ="navajowhite2"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot.allnmds_indsp_temp

nmds.plot.allnmds_indsp_vectors_temp<-nmds.plot.allnmds_indsp_temp+ geom_segment(data = allnmds_indsp_sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = allnmds_indsp_sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label=Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")
nmds.plot.allnmds_indsp_vectors_temp
###



#other ways of doing it....
data.scores <- as.data.frame(scores(allTMT_nmds))  

data.scores$site <- rownames(data.scores)  
data.scores$grp<-trt

ggplot(data=data.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=trt),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=trt,colour=trt),size=4) + 
  theme_mine()
view(OTUallTMTcosq)


OTUallTMT<-dfRowName(OTUallTMTcosq, name ="UNIT")
view(OTUallTMTcosq)
OTUallTMTcosq_HT= merge(HTcol_df, OTUallTMTcosq, by="UNIT")#merging HT but  now need to merge TMT
view(OTUallTMTcosq_HT)
adon.results<-adonis2(OTUallTMTcosq ~ TMT, method="bray", perm=1000)
adon.results<-adonis2(allTMT_site.scrs ~ TMT, method="bray", perm=1000)

write_csv(OTUallTMTco, "OTUallTMTco.csv")



#importing
HT<-read.csv("HT.csv", header = T, strip.white = T, na.strings = "")
#making HTcol a dataframe
HT<-as.data.frame(HT)
view(HT)

TMT<-read.csv("TMT.csv", header = T, strip.white = T, na.strings = "")
TMT<-as.data.frame(TMT)
view(TMT)

OTUallTMTcosq

all.mds <- metaMDS(OTUallTMTcosq)  #OK

data.scores <- as.data.frame(scores(all.mds)) 
data.scores$grp<-TMT


allTMT_nmds <- metaMDS(OTUallTMTco, distance = "bray", autotransform = F)

require(reshape2)
OTUallTMTco$id <- rownames(OTUallTMTco) 
melt(OTUallTMTco)
view(OTUallTMTco)

as.numeric(OTUallTMTco)
class(OTUallTMTco)


all.mds <- metaMDS(all.sites)
data.scores <- as.data.frame(scores(all.mds)) 
data.scores$site <- rownames(data.scores) 
data.scores$grp<-trt
view(trt

as.numeric(unlist(OTUallTMTco))     
view(OTUallTMTco)     
str(OTUallTMTco)


view(all.mds)                                #only defined on a data frame with all numeric-alike variables
view(OTUallTMTco)
class(all.mds)


### kalyn
head(OTUallTMTcosq)
test<-sort(OTUallTMTcosq)
str(OTUallTMTcosq)
df<-OTUallTMTcosq
new_df <- df[ order(row.names(df)), ]
head(new_df)
view(new_df)
test<-df %>% arrange(row.names)
df$name<-row.names
write.csv(OTUallTMTcosq, "OTUallTMTcosq.csv",row.names=T) #this exported the rownames but write_csv didnt
#opened this in xcel and ordered manually.Now bringing back in R and named is OTUallTMTcosqRNO.csv
OTUallTMTcosqRNO<-read.csv("OTUallTMTcosqRNO.csv", header = T, strip.white = T, na.strings = "")
head(OTUallTMTcosqRNO)
view(OTUallTMTcosqRNO)
OTUallTMTcosqRNO$X<-NULL

all.mds <- metaMDS(OTUallTMTcosqRNO)  #OK
data.scores <- as.data.frame(scores(all.mds)) 
data.scores$grp<-TMT
scores(all.mds)
head(all.mds)

adon.result<-adonis2(formula = OTUallTMTcosqRNO ~ TMT, data = TMT, permutations = 99) 
print(adon.result)

adonis2(formula = BC.nmds_allTMT_dist ~ TMT2, data=metaallTMTo, permutations = 999, method = "bray")

BC.nmds_allTMT_dist

#changing pco2 for betadisper
#mutating to include pco2
metaallTMTo_pco2 <- metaallTMTo %>% 
  mutate(TMT2 = TMT)
names(metaallTMTo_pco2)
colnames(metaallTMTo_pco2)[3] <- "pco2" #change colnames of   TMT2 to temp
names(metaallTMTo_pco2)
metaallTMTo_pco2
metaallTMTo_pco2$pco2 <- revalue(metaallTMTo_pco2$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
names(metaallTMTo_pco2)
#mutating to iclude temp
metaallTMTo_pco2 <- metaallTMTo_pco2 %>% 
  mutate(TMT2 = TMT)
names(metaallTMTo_pco2)

#Transform the data sqrt
OTUallTMTcosq<-sqrt(OTUallTMTco)
View(OTUallTMTcosq)
names(OTUallTMTcosq)
#for pairwise analysis cbindg OTUallTMTcosq and 
allTMTpairwise_pco2<-cbind(OTUallTMTcosq, metaallTMTo_pco2)
names(allTMTpairwise_pco2)
allTMTpairwise_pco2


#changing temp for betadisper
#mutating to iclude temp
metaallTMTo_temp <- metaallTMTo %>% 
  mutate(TMT2 = TMT)
names(metaallTMTo_temp)
colnames(metaallTMTo_temp)[3] <- "temp" #change colnames of   TMT2 to temp
names(metaallTMTo_temp)
metaallTMTo_temp
metaallTMTo_temp$temp <- revalue(metaallTMTo_temp$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(metaallTMTo_temp)
#for relative abundance stack plot
allTMTporrelabu<-cbind(OTUallTMTco, metaallTMTo_temp)
head(allTMTporrelabu)
#plot relative abundance
relabunplot<-ggplot(allTMTporrelabu, aes(x = TMT, y = Abundance, fill = TaxallTMT.clean$Class)) + 
  geom_bar(stat = "identity")
relabunplot
#Transform the data sqrt
OTUallTMTcosq<-sqrt(OTUallTMTco)
View(OTUallTMTcosq)
names(OTUallTMTcosq)
#for pairwise analysis cbindg OTUallTMTcosq and 
allTMTpairwise_temp<-cbind(OTUallTMTcosq, metaallTMTo_temp)
names(allTMTpairwise_temp)

