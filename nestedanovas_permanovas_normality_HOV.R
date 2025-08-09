#Demospongiae
#Demospongiae abundance
#Running anova nested within header tanks for abundance of Demospongiae

#abundance

view (demalltmtmergeddfabun)
#edit Unit to UNIT
demalltmtmergeddfabun
colnames(demalltmtmergeddfabun)[2] <- "UNIT" #change colnames of   TMT2 to temp
demalltmtmergeddfabun_HT= merge(HTcol_df, demalltmtmergeddfabun, by="UNIT")
view(demalltmtmergeddfabun_HT)

#adding pco2 and temp as columns
demalltmtmergeddfabun_HT <- demalltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
view(demalltmtmergeddfabun_HT)
colnames(demalltmtmergeddfabun_HT)[6] <- "temp" #change colnames of   TMT2 to temp
names(demalltmtmergeddfabun_HT)
view(demalltmtmergeddfabun_HT)
demalltmtmergeddfabun_HT$temp <- revalue(demalltmtmergeddfabun_HT$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(demalltmtmergeddfabun_HT)

demalltmtmergeddfabun_HT <- demalltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT)
names(demalltmtmergeddfabun_HT)
colnames(demalltmtmergeddfabun_HT)[7] <- "pco2" #change colnames of   TMT2 to temp
names(demalltmtmergeddfabun_HT)
demalltmtmergeddfabun_HT$pco2 <- revalue(demalltmtmergeddfabun_HT$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
view(demalltmtmergeddfabun_HT)

#by treatment
demalltmtmergeddfabun_HT.TMT.base1<-lm(total~TMT/factor(Header), data=demalltmtmergeddfabun_HT) #anova test
demalltmtmergeddfabun_HT.TMT.base1.anovares<-anova(demalltmtmergeddfabun_HT.TMT.base1)
demalltmtmergeddfabun_HT.TMT.base1.anovares

#factorial
demalltmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=demalltmtmergeddfabun_HT) #anova test
demalltmtmergeddfabun_HT.factorial.base1.anovares<-anova(demalltmtmergeddfabun_HT.factorial.base1)
demalltmtmergeddfabun_HT.factorial.base1.anovares

demalltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=demalltmtmergeddfabun_HT) #anova testdemalltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=demalltmtmergeddfabun_HT) #anova test
demalltmtmergeddfabun_HT.factorial.base1.anovares<-anova(demalltmtmergeddfabun_HT.factorial.base1)

hist(residuals(demalltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(demalltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 

#Anderson Darlington normality
ad.test(residuals(demalltmtmergeddfabun_HT.factorial.base1)) 

Anderson-Darling normality test

data:  residuals(demalltmtmergeddfabun_HT.factorial.base1)
A = 0.1591, p-value = 0.942

with(demalltmtmergeddfabun_HT,leveneTest(total ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.3832 0.8988
16  


#export table...
demalltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(demalltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
demalltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(demalltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = demalltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "demalltmtmergeddfabun_HT")
write.csv(demalltmtmergeddfabun_HT.factorial.base1.anovares_df, "demalltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)


#Porifera abundance
#Running anova nested within header tanks for abundance of Keratosa

#abundance

view (poriferaalltmtmergeddfabun)
#edit Unit to UNIT
poriferaalltmtmergeddfabun
colnames(poriferaalltmtmergeddfabun)[2] <- "UNIT" #change colnames of   TMT2 to temp
poriferaalltmtmergeddfabun_HT= merge(HTcol_df, poriferaalltmtmergeddfabun, by="UNIT")
view(poriferaalltmtmergeddfabun_HT)
#by treatment
poralltmtmergeddfabun_HT.TMT.base1<-lm(total~TMT/factor(Header), data=poriferaalltmtmergeddfabun_HT) #anova test
poralltmtmergeddfabun_HT.TMT.base1.anovares<-anova(poralltmtmergeddfabun_HT.TMT.base1)
poralltmtmergeddfabun_HT.TMT.base1.anovares


poralltmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=poriferaalltmtmergeddfabun_HT) #anova test
poralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(poralltmtmergeddfabun_HT.factorial.base1)
poralltmtmergeddfabun_HT.factorial.base1.anovares

poralltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=poralltmtmergeddfabun_HT) #anova testporalltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=poralltmtmergeddfabun_HT) #anova test
poralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(poralltmtmergeddfabun_HT.factorial.base1)

hist(residuals(poralltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(poralltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 
#Anderson Darlington normality
Anderson-Darling normality test

data:  residuals(poralltmtmergeddfabun_HT.factorial.base1)
A = 0.24048, p-value = 0.7481

with(poriferaalltmtmergeddfabun_HT,leveneTest(total ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.4076 0.8838
16 


#export table...
poralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(poralltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
poralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(poralltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = poralltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "poralltmtmergeddfabun_HT")
write.csv(poralltmtmergeddfabun_HT.factorial.base1.anovares_df, "poralltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)


#Keratosa
#Running anova nested within header tanks for abundance of Keratosa

#abundance
keralltmtmerged_abun<-rbind(kerHTHCO2_abund_mergedunit_asdf_total_clean, kerATHCO2_abund_mergedunit_asdf_total_clean, kerAmbient_abund_mergedunit_asdf_total_clean, kerHTACO2_abund_mergedunit_asdf_total_clean)#rbind changes the string of the dataframe. So going to import as .csv and import back in R. This is the only way I am able to do it. I am sure there is a better way. 
keralltmtmerged_abun_df<-as.data.frame(keralltmtmerged_abun) #Had to change from matrix to dataframe
view(keralltmtmerged_abun_df)
write_csv(keralltmtmerged_abun_df, "keralltmtmerged_abun_df.csv")
keralltmtmergeddfabun<-read.csv("keralltmtmerged_abun_df.csv") #import .csv

keralltmtmergeddfabun<-cbind(keralltmtmergeddfabun, Group='Keratosa')
names(keralltmtmergeddfabun)
view(keralltmtmergeddfabun)

#adding pco2 and temp as columns
keralltmtmergeddfabun <- keralltmtmergeddfabun %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
view(keralltmtmergeddfabun)
colnames(keralltmtmergeddfabun)[5] <- "temp" #change colnames of   TMT2 to temp
names(keralltmtmergeddfabun)
view(keralltmtmergeddfabun)
keralltmtmergeddfabun$temp <- revalue(keralltmtmergeddfabun$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(keralltmtmergeddfabun)

keralltmtmergeddfabun <- keralltmtmergeddfabun %>% 
  mutate(TMT2 = TMT)
names(keralltmtmergeddfabun)
colnames(keralltmtmergeddfabun)[6] <- "pco2" #change colnames of   TMT2 to temp
names(keralltmtmergeddfabun)
keralltmtmergeddfabun$pco2 <- revalue(keralltmtmergeddfabun$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
view(keralltmtmergeddfabun)

view(keralltmtmergeddfabun)
#change col names to be the same between dataframes
colnames(keralltmtmergeddfabun)<-c('total','UNIT','TMT', 'Group', 'temp', 'pco2')
view(keralltmtmergeddfabun)
#merging "hapralltmtmergeddfabun" and"HTcol_df based on UNIT

keralltmtmergeddfabun_HT= merge(HTcol_df, keralltmtmergeddfabun, by="UNIT")
view(keralltmtmergeddfabun_HT)
keralltmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=keralltmtmergeddfabun_HT) #anova test
keralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(keralltmtmergeddfabun_HT.factorial.base1)
keralltmtmergeddfabun_HT.factorial.base1.anovares

keralltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=keralltmtmergeddfabun_HT) #anova testkeralltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=keralltmtmergeddfabun_HT) #anova test
keralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(keralltmtmergeddfabun_HT.factorial.base1)

hist(residuals(keralltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(keralltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 
#Anderson Darlington normality
ad.test(residuals(keralltmtmergeddfabun_HT.factorial.base1)) 
Anderson-Darling normality test

data:  residuals(keralltmtmergeddfabun_HT.factorial.base1)
A = 0.85985, p-value = 0.02298
with(keralltmtmergeddfabun_HT,leveneTest(logtranadd1 ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.8998   0.53
16 

#levene's test
with(keralltmtmergeddfabun_HT,leveneTest(total ~ TMT))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  3  1.4783 0.2507
20     


#export table...
keralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(keralltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
keralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(keralltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = keralltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "keralltmtmergeddfabun_HT")
write.csv(keralltmtmergeddfabun_HT.factorial.base1.anovares_df, "keralltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)


#Running anova nested within header tanks for abundance of Homoscleromorpha
view(homoalltmtmergeddfabun)
#change col names to be the same between dataframes
colnames(homoalltmtmergeddfabun)<-c('total','UNIT','TMT', 'Group', 'temp', 'pco2')
view(homoalltmtmergeddfabun)
#merging "hapralltmtmergeddfabun" and"HTcol_df based on UNIT

homoalltmtmergeddfabun_HT= merge(HTcol_df, homoalltmtmergeddfabun, by="UNIT")
view(homoalltmtmergeddfabun_HT)
homoalltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=homoalltmtmergeddfabun_HT) #anova test

anova(homoalltmtmergeddfabun_HT.base1)

Analysis of Variance Table

Response: total
Df Sum Sq Mean Sq F value  Pr(>F)  
TMT                 3 2689.1  896.37  3.7213 0.03333 *
  TMT:factor(Header)  4  434.8  108.71  0.4513 0.77006  
Residuals          16 3854.0  240.87                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
view(HTcol_df)

#homoalltmtmergeddfabun_HT.base1
hist(residuals(homoalltmtmergeddfabun_HT.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(homoalltmtmergeddfabun_HT.base1) #qqresiduals shows variance to be somewhat stable. 
#Anderson Darlington normality
ad.test(residuals(homoalltmtmergeddfabun_HT.base1)) 
Anderson-Darling normality test
#results show that TMT had a significant impact on abundanct of homo but not header tanks. 

homoalltmtmergeddfabun_HT$logtranadd1<-log (homoalltmtmergeddfabun_HT$total+1) #log transformation
homoalltmtmergeddfabun_HT.base1<-lm(logtranadd1~TMT/factor(Header), data=homoalltmtmergeddfabun_HT) #anova test
anova(homoalltmtmergeddfabun_HT.base1)
Analysis of Variance Table

Response: logtranadd1
Df Sum Sq Mean Sq F value   Pr(>F)   
TMT                 3 36.004 12.0014  7.8408 0.001926 **
  TMT:factor(Header)  4  0.958  0.2395  0.1564 0.957238   
Residuals          16 24.490  1.5306                    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1      
      

homoalltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=homoalltmtmergeddfabun_HT) #anova test
homoalltmtmergeddfabun_HT.factorial.base1.anovares<-anova(homoalltmtmergeddfabun_HT.factorial.base1)
Analysis of Variance Table

Response: logtranadd1
Df  Sum Sq Mean Sq F value   Pr(>F)    
pco2                      1 31.2867 31.2867 20.4404 0.000348 ***
  temp                      1  3.7523  3.7523  2.4515 0.136978    
pco2:temp                 1  0.9653  0.9653  0.6307 0.438733    
pco2:temp:factor(Header)  4  0.9578  0.2395  0.1564 0.957238    
Residuals                16 24.4901  1.5306                     
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#export table...

homoalltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(homoalltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df

#then add the names and data or variables...
homoalltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(homoalltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = homoalltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "homoalltmtmergeddfabun_HT")
write.csv(homoalltmtmergeddfabun_HT.factorial.base1.anovares_df, "homoalltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)

hist(residuals(homoalltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(homoalltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 
#Anderson Darlington normality
ad.test(residuals(homoalltmtmergeddfabun_HT.factorial.base1)) 
Anderson-Darling normality test

data:  residuals(homoalltmtmergeddfabun_HT.factorial.base1)
A = 0.85985, p-value = 0.02298
with(homoalltmtmergeddfabun_HT,leveneTest(logtranadd1 ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
      Df F value Pr(>F)
group  7  0.8998   0.53
      16 


      
#parametric with tranformation and do a factorial analysis...
view(homoalltmtmergeddfabun_HT)
homocumabun.ra.logtranadd1_HT<-lme(logtranadd1 ~ pco2*temp/factor(Header), random=~ 1|UNIT, data=homoalltmtmergeddfabun_HT, method = "REML")
homocumabun.ra.logtranadd1.aov<-aov(homocumabun.ra.logtranadd1, data=homoalltmtmergeddfabun_HT)
homocumabun.ra.logtranadd1.aov.sum<-summary(homocumabun.ra.logtranadd1.aov)
summary(homocumabun.ra.logtranadd1.aov)

#Tetractinellida abundance
#Running anova nested within header tanks for abundance of Demospongiae

#abundance

view (tetralltmtmergeddfabun)
#edit Unit to UNIT
tetralltmtmergeddfabun
colnames(tetralltmtmergeddfabun)[2] <- "UNIT" #change colnames of   TMT2 to temp
tetralltmtmergeddfabun_HT= merge(HTcol_df, tetralltmtmergeddfabun, by="UNIT")
view(tetralltmtmergeddfabun_HT)

#adding pco2 and temp as columns
tetralltmtmergeddfabun_HT <- tetralltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
view(tetralltmtmergeddfabun_HT)
colnames(tetralltmtmergeddfabun_HT)[5] <- "temp" #change colnames of   TMT2 to temp
names(tetralltmtmergeddfabun_HT)
view(tetralltmtmergeddfabun_HT)
tetralltmtmergeddfabun_HT$temp <- revalue(tetralltmtmergeddfabun_HT$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(tetralltmtmergeddfabun_HT)

tetralltmtmergeddfabun_HT <- tetralltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT)
names(tetralltmtmergeddfabun_HT)
colnames(tetralltmtmergeddfabun_HT)[6] <- "pco2" #change colnames of   TMT2 to temp
names(tetralltmtmergeddfabun_HT)
tetralltmtmergeddfabun_HT$pco2 <- revalue(tetralltmtmergeddfabun_HT$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
view(tetralltmtmergeddfabun_HT)

#by treatment
tetralltmtmergeddfabun_HT.TMT.base1<-lm(total~TMT/factor(Header), data=tetralltmtmergeddfabun_HT) #anova test
tetralltmtmergeddfabun_HT.TMT.base1.anovares<-anova(tetralltmtmergeddfabun_HT.TMT.base1)
tetralltmtmergeddfabun_HT.TMT.base1.anovares

#factorial
tetralltmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=tetralltmtmergeddfabun_HT) #anova test
tetralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(tetralltmtmergeddfabun_HT.factorial.base1)
tetralltmtmergeddfabun_HT.factorial.base1.anovares

tetralltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=tetralltmtmergeddfabun_HT) #anova testtetralltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=tetralltmtmergeddfabun_HT) #anova test
tetralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(tetralltmtmergeddfabun_HT.factorial.base1)

hist(residuals(tetralltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(tetralltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 

#Anderson Darlington normality
ad.test(residuals(tetralltmtmergeddfabun_HT.factorial.base1)) 

Anderson-Darling normality test

data:  residuals(tetralltmtmergeddfabun_HT.factorial.base1)
A = 0.2524, p-value = 0.7072

with(tetralltmtmergeddfabun_HT,leveneTest(total ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.9156 0.5195
16    


#export table...
tetralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(tetralltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
tetralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(tetralltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = tetralltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "tetralltmtmergeddfabun_HT")
write.csv(tetralltmtmergeddfabun_HT.factorial.base1.anovares_df, "tetralltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)


#Poecilosclerida abundance
#Running anova nested within header tanks for abundance of Demospongiae

#abundance poecralltmtmergeddfabun

view (poecralltmtmergeddfabun)
#edit Unit to UNIT
poecralltmtmergeddfabun
colnames(poecralltmtmergeddfabun)[2] <- "UNIT" #change colnames of   TMT2 to temp
poecralltmtmergeddfabun_HT= merge(HTcol_df, poecralltmtmergeddfabun, by="UNIT")
view(poecralltmtmergeddfabun_HT)

#adding pco2 and temp as columns
poecralltmtmergeddfabun_HT <- poecralltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
view(poecralltmtmergeddfabun_HT)
colnames(poecralltmtmergeddfabun_HT)[5] <- "temp" #change colnames of   TMT2 to temp
names(poecralltmtmergeddfabun_HT)
view(poecralltmtmergeddfabun_HT)
poecralltmtmergeddfabun_HT$temp <- revalue(poecralltmtmergeddfabun_HT$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(poecralltmtmergeddfabun_HT)

poecralltmtmergeddfabun_HT <- poecralltmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT)
names(poecralltmtmergeddfabun_HT)
colnames(poecralltmtmergeddfabun_HT)[6] <- "pco2" #change colnames of   TMT2 to temp
names(poecralltmtmergeddfabun_HT)
poecralltmtmergeddfabun_HT$pco2 <- revalue(poecralltmtmergeddfabun_HT$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
view(poecralltmtmergeddfabun_HT)

#by treatment
poecralltmtmergeddfabun_HT.TMT.base1<-lm(total~TMT/factor(Header), data=poecralltmtmergeddfabun_HT) #anova test
poecralltmtmergeddfabun_HT.TMT.base1.anovares<-anova(poecralltmtmergeddfabun_HT.TMT.base1)
poecralltmtmergeddfabun_HT.TMT.base1.anovares

#factorial
poecralltmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=poecralltmtmergeddfabun_HT) #anova test
poecralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(poecralltmtmergeddfabun_HT.factorial.base1)
poecralltmtmergeddfabun_HT.factorial.base1.anovares

poecralltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=poecralltmtmergeddfabun_HT) #anova testpoecralltmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=poecralltmtmergeddfabun_HT) #anova test
poecralltmtmergeddfabun_HT.factorial.base1.anovares<-anova(poecralltmtmergeddfabun_HT.factorial.base1)

hist(residuals(poecralltmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(poecralltmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 

#Anderson Darlington normality
ad.test(residuals(poecralltmtmergeddfabun_HT.factorial.base1)) 

Anderson-Darling normality test

data:  residuals(poecralltmtmergeddfabun_HT.factorial.base1)
A = 0.50509, p-value = 0.1833


with(poecralltmtmergeddfabun_HT,leveneTest(total ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.1934 0.9826
16   


#export table...
poecralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(poecralltmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
poecralltmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(poecralltmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = poecralltmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "poecralltmtmergeddfabun_HT")
write.csv(poecralltmtmergeddfabun_HT.factorial.base1.anovares_df, "poecralltmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)

#Suberitida abundance
#Running anova nested within header tanks for abundance of Demospongiae

#abundance suballtmtmergeddfabun

view (suballtmtmergeddfabun)
#edit Unit to UNIT
suballtmtmergeddfabun
colnames(suballtmtmergeddfabun)[2] <- "UNIT" #change colnames of   TMT2 to temp
suballtmtmergeddfabun_HT= merge(HTcol_df, suballtmtmergeddfabun, by="UNIT")
view(suballtmtmergeddfabun_HT)

#adding pco2 and temp as columns
suballtmtmergeddfabun_HT <- suballtmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT) #copying the same column twice to make a factorial design
view(suballtmtmergeddfabun_HT)
colnames(suballtmtmergeddfabun_HT)[5] <- "temp" #change colnames of   TMT2 to temp
names(suballtmtmergeddfabun_HT)
view(suballtmtmergeddfabun_HT)
suballtmtmergeddfabun_HT$temp <- revalue(suballtmtmergeddfabun_HT$temp, c("ATACO2"="Present","ATHCO2"="Present", "HTACO2"="Future", "HTHCO2"="Future"))
names(suballtmtmergeddfabun_HT)

suballtmtmergeddfabun_HT <- suballtmtmergeddfabun_HT %>% 
  mutate(TMT2 = TMT)
names(suballtmtmergeddfabun_HT)
colnames(suballtmtmergeddfabun_HT)[6] <- "pco2" #change colnames of   TMT2 to temp
names(suballtmtmergeddfabun_HT)
suballtmtmergeddfabun_HT$pco2 <- revalue(suballtmtmergeddfabun_HT$pco2, c("ATACO2"="Present","ATHCO2"="Future", "HTACO2"="Present", "HTHCO2"="Future"))
view(suballtmtmergeddfabun_HT)

#by treatment
suballtmtmergeddfabun_HT.TMT.base1<-lm(total~TMT/factor(Header), data=suballtmtmergeddfabun_HT) #anova test
suballtmtmergeddfabun_HT.TMT.base1.anovares<-anova(suballtmtmergeddfabun_HT.TMT.base1)
suballtmtmergeddfabun_HT.TMT.base1.anovares

#factorial
suballtmtmergeddfabun_HT.factorial.base1<-lm(total~pco2*temp/factor(Header), data=suballtmtmergeddfabun_HT) #anova test
suballtmtmergeddfabun_HT.factorial.base1.anovares<-anova(suballtmtmergeddfabun_HT.factorial.base1)
suballtmtmergeddfabun_HT.factorial.base1.anovares

suballtmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=suballtmtmergeddfabun_HT) #anova testsuballtmtmergeddfabun_HT.factorial.base1<-lm(logtranadd1~pco2*temp/factor(Header), data=suballtmtmergeddfabun_HT) #anova test
suballtmtmergeddfabun_HT.factorial.base1.anovares<-anova(suballtmtmergeddfabun_HT.factorial.base1)

hist(residuals(suballtmtmergeddfabun_HT.factorial.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(suballtmtmergeddfabun_HT.factorial.base1) #qqresiduals shows variance to be somewhat stable. 

#Anderson Darlington normality
ad.test(residuals(suballtmtmergeddfabun_HT.factorial.base1)) 

Anderson-Darling normality test

data:  residuals(suballtmtmergeddfabun_HT.factorial.base1)
A = 0.12157, p-value = 0.9858


with(suballtmtmergeddfabun_HT,leveneTest(total ~ TMT/factor(Header)))
#Levene's Test for Homogeneity of Variance (center = median)
Df F value Pr(>F)
group  7  0.4847  0.832
16    


#export table...
suballtmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(suballtmtmergeddfabun_HT.factorial.base1.anovares)#convert to df
#then add the names and data or variables...
suballtmtmergeddfabun_HT.factorial.base1.anovares_df <- tidy(suballtmtmergeddfabun_HT.factorial.base1.anovares) %>% 
  mutate(data.name = suballtmtmergeddfabun_HT.factorial.base1.anovares$data.name, 
         dataset = "suballtmtmergeddfabun_HT")
write.csv(suballtmtmergeddfabun_HT.factorial.base1.anovares_df, "suballtmtmergeddfabun_HT.factorial.base1.anovares_df.csv", row.names = F)






#Running anova nested within header tanks for abundance of Haplosclerida
view(HTcol_df)
#merging df based on UNIT
colnames(hapralltmtmergeddfabun)<-c('total','UNIT','TMT')
hapralltmtmergeddfabun_HT= merge(HTcol_df, hapralltmtmergeddfabun, by="UNIT")
hapralltmtmergeddfabun_HT.base1<-lm(total~TMT/factor(Header), data=hapralltmtmergeddfabun_HT) #anova test
anova(hapralltmtmergeddfabun_HT.base1)
view(hapralltmtmergeddfabun)
view(HTcol_df)

Response: total
Df Sum Sq Mean Sq F value  Pr(>F)  
TMT                 3   3125 1041.67  3.0990 0.05643 .
TMT:factor(Header)  4    615  153.75  0.4574 0.76579  
Residuals          16   5378  336.13                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#running nested analysis of variance 
#merging df based on UNIT
poralltmtmergeddfrichness_HT= merge(HTcol_df, poralltmtmergeddfrichness, by="UNIT")
view(poralltmtmergeddfrichness_HT)
view(metaallTMTo)

Analysis of Variance Table

Response: Observed
Df  Sum Sq Mean Sq F value Pr(>F)
TMT                 3  15.167  5.0556  0.7490 0.5387
TMT:factor(Header)  4  54.667 13.6667  2.0247 0.1393
Residuals          16 108.000  6.7500   


#nested anova on diversity for all sponges 
view(poralltmtmergeddfrichness_HT)
poralltmtmergeddfrichness_HT.base1<-lm(Observed~TMT/factor(Header), data=poralltmtmergeddfrichness_HT)
anova(poralltmtmergeddfrichness_HT.base1)
Analysis of Variance Table

Analysis of Variance Table

Response: Observed
Df  Sum Sq Mean Sq F value Pr(>F)
TMT                 3  15.167  5.0556  0.7490 0.5387
TMT:factor(Header)  4  54.667 13.6667  2.0247 0.1393
Residuals          16 108.000  6.7500    
poralltmtmergeddfrichness.base1.aov<-anova(poralltmtmergeddfrichness.base1)


hist(residuals(poralltmtmergeddfrichness.base1))
layout(matrix(c(1,2,3,4),2,2))
plot(poralltmtmergeddfrichness.base1) #qqresiduals shows variance to be somewhat stable. 
#Anderson Darlington normality
ad.test(residuals(poralltmtmergeddfrichness.base1))
Anderson-Darling normality test

data:  residuals(poralltmtmergeddfrichness.base1)
A = 0.2615, p-value = 0.6755

#HOV assumption-Levene's test on TMT and it passes homogeneity of variance
with(poralltmtmergeddfrichness,leveneTest(Observed ~ TMT))
#Levene's Test for Homogeneity of Variance (center = median)
#Df F value Pr(>F)
#group  3  1.9833 0.1489
#20    





#Adonis nested (not sure we can do this)
BC.nmds_allTMT_dist_adonis<-adonis2(formula = BC.nmds_allTMT_dist ~ TMT/factor(Header), data = metaallTMTo_HT, permutations = 1000)
BC.nmds_allTMT_dist_adonis
Df SumOfSqs      R2      F   Pr(>F)    
Model     3  0.60906 0.23195 2.0133 0.000999 ***
  Residual 20  2.01677 0.76805                    
Total    23  2.62583 1.00000                    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pairwise comparisons between Management levels
metaallTMTo$TMT=as.factor(metaallTMTo$TMT)
# Pairwise comparisons between Management levels
pwnmdsTMT <- adonis_pairwise(x = metaallTMTo, dd = vegdist(BC.nmds_allTMT_dist), group.var = "TMT")
pwnmdsTMT
class(pwnmdsTMT)