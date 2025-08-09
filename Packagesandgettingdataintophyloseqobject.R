
```{r setup chunk, setup, include = FALSE, cache=FALSE, message=FALSE, warning=FALSE}
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('xtable')) install.packages('xtable'); library('RColorBrewer')
if (!require('utils')) install.packages('utils'); library('utils')
if (!require('nlme')) install.packages('nlme'); library('nlme')
if (!require('stats')) install.packages('stats'); library('stats')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('devtools')) install.packages('devtools'); library('devtools')
if (!require('multcomp')) install.packages('multcomp'); library('multcomp')
if (!require('multcompView')) install.packages('multcompView'); library('multcompView')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('lmerTest')) install.packages('lmerTest'); library('lmerTest')
if (!require('lsmeans')) install.packages('lsmeans'); library('lsmeans')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('knitr')) install.packages('knitr'); library('knitr')
#download the following packages
if (!require('picante')) install.packages('picante'); library('picante')
if (!require('broom')) install.packages('broom'); library('broom')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('vegan')) install.packages('vegan'); library('vegan')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('Hmisc')) install.packages('Hmisc'); library('Hmisc')
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('viridis')) install.packages('viridis'); library('viridis')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('plotrix')) install.packages('plotrix'); library('plotrix')
if (!require('Rmisc')) install.packages('Rmisc'); library('Rmisc')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('caroline')) install.packages('caroline'); library('caroline')
if (!require('gganimate')) install.packages('gganimate'); library('gganimate')
if (!require('gifski')) install.packages('gifski'); library('gifski')
if (!require('png')) install.packages('png'); library('png')
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('nortest')) install.packages('nortest'); library('nortest')
if (!require('car')) install.packages('car'); library('car')
if (!require('heatmaply')) install.packages('heatmaply'); library('heatmaply')
if (!require('reshape')) install.packages('reshape'); library('reshape')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('textshape')) install.packages('textshape'); library('textshape')
if (!require('phyloseq')) install.packages('phyloseq'); library('phyloseq')
if (!require('utils')) install.packages('utils'); library('utils')
if (!require('tibble')) install.packages('tibble'); library('tibble')
if (!require('ggConvexHull')) install.packages('ggConvexHull'); library('ggConvexHull')
if (!require('remotes')) install.packages('remotes'); library('remotes')
if (!require('RVAideMemoire')) install.packages('RVAideMemoire'); library('RVAideMemoire')
if (!require('pscl')) install.packages('pscl'); library('pscl')
if (!require('performance')) install.packages('performance'); library('performance')
if (!require('networkD3')) install.packages('networkD3'); library('networkD3')
if (!require('svglite')) install.packages('svglite'); library('svglite')
if (!require('webshot')) install.packages('webshot'); library('webshot')
if (!require('speedyseq')) install.packages('speedyseq'); library('speedyseq')
install.packages("remotes")
library(remotes)
install_version("webshot", "0.5.4")
if (!require('webshot2')) install.packages('webshot2'); library('webshot2')
install.packages("webshot", version='0.5.4');library('webshot')
webshot::install_phantomjs(force=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install_github("Russel88/MicEco", force=TRUE)
install_github("adrientaudiere/Miscmetabar",force=TRUE)
library('MiscMetabar')
library('MicEco')

```


```{r setup, include=FALSE}
#packages installed through BiocManager or devetools or remotes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rhdf5lib")
devtools::install_github("cmartin/ggConvexHull", force = TRUE)
BiocManager::install("limma")
BiocManager::install("phyloseq", force = TRUE)
devtools::install_github("dkahle/ggmap")
devtools::install_github("mikemc/speedyseq", force = TRUE)#helps edit tax table in phyloseq object
devtools::install_github("adrientaudiere/Miscmetabar", force = TRUE)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true") # need to set this to true in order for package to install
remotes::install_github("Russel88/MicEco") #Produces very nice venn diagrams. Ran into and issue with this and kept on getting the error message. sing GitHub PAT from the git credential store.Error: Failed to install 'unknown package' from GitHub:HTTP error 401.Bad credentials. was able to solve it by using "gitcreds::gitcreds_delete()" and then presseing "2. Delte these credentials".
devtools::install_github("vmikk/metagMisc") #install mEtagMisx allows to export phyloseq objects to df
devtools::install_github("tomwenseleers/export")


```

#Functions
```{r setup, include=FALSE}
#Add function metadata
add_metadata <- function(x, metad, xid = NULL, mid = NULL, drop_mid = T){
  
  ## Data validation
  if(is.null(xid) | is.null(mid)){ stop("Error: row identifiers should be provided.\n") }
  if(!length(xid) == 1){ stop("Error: row identifiers should be provided as a single character string, e.g. 'SampleID'.\n") }
  if(!length(mid) == 1){ stop("Error: row identifiers should be provided as a single character string, e.g. 'SampleID'.\n") }
  
  ## Extract metadata from phyloseq object
  if(class(metad) %in% "phyloseq"){
    if(is.null(phyloseq::sample_data(metad, errorIfNULL = F))){
      stop("Error: sample_data is missing from the phyloseq object 'metad'.\n")
    } else {
      metad <- as(object = phyloseq::sample_data(metad), Class = "data.frame")
    }
  }
  
  ## Data validation
  if(!xid %in% colnames(x)){ stop("Error: '", xid, "' column is missing in the main data.\n", sep="") }
  if(!mid %in% colnames(metad)){ stop("Error: '", mid, "' column is missing in the metadata.\n", sep="") }
  
  if(nrow(x) != length(unique(x[, xid]))){ stop("Error: Row identifiers are not unique in 'x'.\n") }
  if(nrow(metad) != length(unique(metad[, mid]))){ stop("Error: Row identifiers are not unique in 'metad'.\n") }
  
  
  ## Match metadata to the main data
  mm <- match(x = x[, xid], table = metad[, mid])
  
  ## Reorder metadata
  metad <- metad[mm, ]
  
  ## Remove sample ID column from metadata
  if(drop_mid == TRUE){
    metad <- metad[, -which(colnames(metad) == mid)]
  }
  
  ## Merge data with metadata
  res <- cbind(x, metad)
  
  return(res)
}

#add function dfRowName
dfRowName <- function(x, name = "Rows", stringsAsFactors = FALSE){
  res <- data.frame(rownames(x), x, stringsAsFactors = stringsAsFactors)
  colnames(res)[1] <- name
  rownames(res) <- NULL
  return(res)
}
```

## Importing and formatting data xcel sheets to R 
#the number "0" for OTUs 9, 19, and 23 appear on the left and not the right which seems to causing a problem when running "data.matrix" command because it ends up adding a number "2" instead of "0"
#here is the problem!!!! as.matrix should be data.matrix. these columns are considered "chr" and not "int" so need to change them.by using integer(). Best way to deal with thi
#is by using header = T, strip.white = T, na.strings = "") in the import step!!!
#porifera<-read.csv("OTUtable.csv", header = T, strip.white = T, na.strings = "")
```{r setup, include=FALSE}
porifera<-read.csv("OTUtable.csv", header = T, strip.white = T, na.strings = "")#better to import as .csv this creates some to be considered characters. so we need to change to integer.
str(porifera)#import as .csv this creates some to be considered characters. so we need to change to integer.
View(porifera) 
porifera$OTU9 <- as.integer(porifera$OTU9)
porifera$OTU19 <- as.integer(porifera$OTU19)
porifera$OTU23 <- as.integer(porifera$OTU23)
porifera$OTU85 <- as.integer(porifera$OTU85)
porifera$OTU86 <- as.integer(porifera$OTU86)
porifera=porifera[,!(names(porifera) %in% c("X"))] #deleting column "X" because its useless
View(porifera)
porifera_df<-as.data.frame(porifera) # Making a dataframe to move forward
str(porifera_df)
porifera_df$Recovery<-as.factor(porifera_df$Recovery) # This is a dummy factor and thus need to remove it from being read as a number
# # changing NA to 0

## But first need to make sure all columns that might be considered numeric are not
porifera_df_numeric<-porifera_df
str(porifera_df_numeric)
names(porifera_df_numeric)


## Changing Date, Unit, Plate, Side, Recovery, SampleTime, PlateTrue to factors
metaNeeds<-c("Date","TMT","UNIT", "PLATE", "SIDE", "Recovery")
porifera_df_numeric[,c(metaNeeds)]<-sapply(porifera_df_numeric[,c(metaNeeds)], as.factor)
porifera_df[is.na(porifera_df)] <- 0
names(porifera_df)

porifera_df$Date<-as.factor(porifera_df$Date)
porifera_df$TMT<-as.factor(porifera_df$TMT)
porifera_df$UNIT<-as.factor(porifera_df$UNIT)
porifera_df$PLATE<-as.factor(porifera_df$PLATE)
porifera_df$SIDE<-as.factor(porifera_df$SIDE)
porifera_df$Recovery<-as.factor(porifera_df$Recovery)                          


# Adding time as 'sample'
porifera_df$SampleTime <- paste(porifera_df$Date,porifera_df$TMT,porifera_df$UNIT,porifera_df$PLATE,porifera_df$SIDE,sep="_")
names(porifera_df)
porifera_df$SampleTime<-as.factor(porifera_df$SampleTime)

#AddingUnit to time
porifera_df$UNITTime<-paste(porifera_df$Date,porifera_df$UNIT,sep="_")
names(porifera_df)
porifera_df$UNITTime<-as.factor(porifera_df$UNITTime)

#Adding TMT to Time
porifera_df$TMTTime<-paste(porifera_df$Date,porifera_df$TMT,sep="_")
names(porifera_df)
porifera_df$TMTTime<-as.factor(porifera_df$TMTTime)

#Adding Treatment to Unit
porifera_df$Unit_Treatment<-paste(porifera_df$TMT, porifera_df$UNIT, sep="")
names(porifera_df)
porifera_df$Unit_Treatment<-as.factor(porifera_df$Unit_Treatment)

# Adding "true" plate column
porifera_df$PlateTrue<-porifera_df$PLATE
porifera_df$PlateTrue[porifera_df$PlateTrue == 1]<-1
porifera_df$PlateTrue[porifera_df$PlateTrue == 2]<-1
porifera_df$PlateTrue[porifera_df$PlateTrue == 3]<-2
porifera_df$PlateTrue[porifera_df$PlateTrue == 4]<-2
porifera_df$PlateTrue[porifera_df$PlateTrue == 5]<-3
porifera_df$PlateTrue[porifera_df$PlateTrue == 6]<-3

porifera_df$PlateTrue<-as.factor(porifera_df$PlateTrue)

#Adding side to  plateTrue
porifera_df$PlateTrueSide<-paste(porifera_df$PlateTrue, porifera_df$SIDE, sep="")
names(porifera_df)
porifera_df$PlateTrueSide<-as.factor(porifera_df$PlateTrueSide)

###### Exporting Cleaned porifera_df file #####
# this will save the file to your working directory
row.names(porifera_df) <- porifera_df$SampleTime #make SampleTime rowname for OTU table and metadata file.

#Making otu file with new sample names as row names
porifera_df_OTU<-porifera_df
names(porifera_df_OTU)
porifera_df_OTU<-porifera_df_OTU[,-(1:6)]
names(porifera_df_OTU)
porifera_df_OTU<-porifera_df_OTU[,-(94:99)]
names(porifera_df_OTU)
write.csv(porifera_df_OTU, "porifera_df_OTU.csv", row.names = F)

#making metadata file as well. 
porifera_df_metadata<-porifera_df
names(porifera_df_metadata)
porifera_df_metadata<-porifera_df_metadata[,-(7:99)]
names(porifera_df_metadata)
write.csv(porifera_df_metadata, "porifera_df_metadata.csv", row.names = F)

#making OTUclass file as well. 
OTUClass<-read_excel("OTUhporder.xlsx", 1)
head(OTUClass)
OTUClass_df<-as.data.frame(OTUClass)
row.names(OTUClass_df) <- OTUClass_df$OTU #make SampleTime rowname for OTU table and metadata file.
names(OTUClass_df)
```

#getting data in phyloseq
```{r setup, include=FALSE}
#Preparing files for phyloseq
#Make dataframes matrices 
#making data file a matrix
porifera_df_OTU.mat<-data.matrix(porifera_df_OTU, rownames.force = NA)
row.names(porifera_df_OTU.mat)
class(porifera_df_OTU.mat)
dim(porifera_df_OTU.mat)

#making classification file a matrices
OTUClass_df

OTUClass_df.mat <- as.matrix(OTUClass_df) #was using OTUClass_df.mat <- data.matrix(OTUClass_df) #having issues with this command becuase it converts chr to numbers but the OTUClass_df needs to be in matrix form inorder for the TAX to accept into phyloseq

class(OTUClass_df.mat)
dim(OTUClass_df.mat)
View(OTUClass_df.mat)

otumat = porifera_df_OTU.mat
taxmat = OTUClass_df.mat
head(taxmat)
head(OTUClass_df.mat)

sampledata=porifera_df_metadata
sampledata=sample_data(sampledata)
sampledata

OTU = otu_table(otumat, taxa_are_rows = FALSE)
TAX =tax_table(taxmat)
physeq = phyloseq(OTU, TAX, sampledata)
physeq

ntaxa (physeq)
nsamples (physeq)
rank_names (physeq)

#transform all variables to factors
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)
```