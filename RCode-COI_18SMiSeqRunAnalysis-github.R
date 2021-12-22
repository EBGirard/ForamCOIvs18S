
rm(list=ls())

#package needed
library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(car)
library(phytools)
library(ape)
library(Biostrings)
library(ggtree)
library(treeio)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%


#-----------Rarefaction curves COI --------------------------------

df_pro <- read.csv("~/Downloads/Dataset/ASVtable-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/COI-raw-reads-table.csv", sep = ";")

df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[,c(135,134)]
df_raw1 <- df_raw[c(1,2)]

df_all <- merge(df_raw1, df_pro1, by = "Sample_name")
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100

summary(df_all$ratio)

#retain samples with ratio > 10% 
sample_to_analyze <- levels(droplevels(df_all$Sample_name[df_all$ratio > 10]))

#COI create dataset with different read filter threshold for sampling curve

data <- read.csv("~/Downloads/Dataset/ASVtable-new.csv", sep = ";", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/MetadataCOIsamples20210726.csv", sep = ";")

#remove samples with read ratio < 10%
keep <- sample_to_analyze
data <- data[,(names(data) %in% keep)]
data <- data[rowSums(data) > 0, ]

#create dataset for sampling curve
cutoff0 <- rep(0,nrow(data))
cutoff0.1 <- rep(0,nrow(data))
cutoff0.5 <- rep(0,nrow(data))
cutoff1 <- rep(0,nrow(data))
cutoff1.5 <- rep(0,nrow(data))
cutoff2 <- rep(0,nrow(data))
cutoff2.5 <- rep(0,nrow(data))
cutoff3 <- rep(0,nrow(data))
cutoff5 <- rep(0,nrow(data))
cutoff10 <- rep(0,nrow(data))
df_curve <- data.frame(cutoff0, cutoff0.1,cutoff0.5,cutoff1, cutoff1.5, cutoff2, cutoff2.5, cutoff3,cutoff5, cutoff10)
df_curve$otuid <- rownames(data)
df_curve <- df_curve[, c(ncol(df_curve),1:(ncol(df_curve)-1))]

#_____________________________run loop from here________________________________________

cuts <- c(0, 0, 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1)

for (c in 2:ncol(df_curve)){
  
  otudata_tfilt <- data.frame(t(data))
  otudf <- otudata_tfilt
  otudata_tfilt <- as.data.frame(t(apply(otudata_tfilt, 1, function(x) x/sum(x))))
  rowSums(otudata_tfilt)
  #remove asv < 0.5% of reads of the total reads but keeping reads number
  for (m in 1:ncol(otudf)) {
    
    for (v in 1:nrow(otudf)) {
      
      if (otudata_tfilt[v,m] < cuts[c]) { #to account for cross contamination and barcode switching during sequencing
        
        otudf[v,m] <- 0
      }}}
  rowSums(otudf)
  colSums(otudf)
  otudf <- otudf[,colSums(otudf) > 0]
  #join metadata with asv data
  otudf <- data.frame(otudf)
  otudf$baseclear_id <- rownames(otudf)
  df10 <- merge(metadata[, c(1,2,7,8,10)], otudf, by = "baseclear_id")
  #add island number only
  splitnumber <- str_split_fixed(df10$field.nmbr, "-", 2)
  df10$island <-  splitnumber[,1]
  df10<- df10[,c(1:5,ncol(df10),6:(ncol(df10)-1))]
  
  df <- df10
  df <- melt(df, id.vars = c("baseclear_id","mastermix","registr.code","field.nmbr","taxon","island"), 
             value.name = "reads", variable.name = "Name")
  df <- df[,c(1,2,3,4,5,6,7,8)]
  df <- df[df$reads > 0, ]
  
  cutoffasv <- levels(droplevels(df$Name))
  readsum <- df %>% group_by(Name) %>% summarise(totread = sum(reads))
  
  for (x in 1:length(cutoffasv)){
    
    for(i in 1:nrow(df_curve)) {
      
      if (cutoffasv[x] == df_curve$otuid[i]) {
        
        df_curve[i,c] <- readsum$totread[x]
        
      }
    }
  }
  
}

df55 <- df_curve

rownames(df55) <- df55$otuid
df55 <- t(df55[,-1])


#plot rarefaction curve
S <- specnumber(df55) # observed number of species
raremax <- min(rowSums(df55))
Srare <- rarefy(df55, raremax)

data1 <- data
for (j in 2:ncol(data1)){
  for(i in 1:nrow(data1)) {
    if (data1[i,j] > 1) {
      
      data1[i,j] <- 1
    }
  }
  
}

colSums(data1[2:ncol(data1)])
max(colSums(data[2:ncol(data)]))

par(mfrow = c(1,2))
plot(colSums(data1[2:ncol(data1)]), colSums(data[2:ncol(data)])/1000, xlab = "Number of ASVs per sample",
     ylab = "Number of reads per sample (in thousands)", main = "Before any threshold applied")
rarecurve(df55, step = 20, col = "blue", cex = 0.6, 
          xlab = "Total number of reads", ylab = "Number of ASVs", xlim = c(0, 50000))




#-----------Rarefaction curves 18S --------------------------------
df_pro <- read.csv("~/Downloads/Dataset/ASVtable-18S-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/18S-raw-reads-table.csv", sep = ";")

df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[df_pro_t$processed_read_sum > 1000,c(1385,1384)]
df_raw1 <- df_raw[c(1,2)]



df_all <- merge(df_raw1, df_pro1, by = "Sample_name")
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100

summary(df_all$ratio)

#retain samples with ratio > 10% 
sample_to_analyze <- levels(droplevels(df_all$Sample_name[df_all$ratio > 10]))

#18S create dataset with different read filter threshold for sampling curve

data <- read.csv("~/Downloads/Dataset/ASVtable-18S-new.csv", sep = ";", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/Metadata18Ssamples20210726.csv", sep = ";")

#remove samples with read ratio < 10%
keep <- sample_to_analyze
data <- data[keep]
data <- data[rowSums(data) > 0, ]

#create dataset for sampling curve
cutoff0 <- rep(0,nrow(data))
cutoff0.1 <- rep(0,nrow(data))
cutoff0.5 <- rep(0,nrow(data))
cutoff1 <- rep(0,nrow(data))
cutoff1.5 <- rep(0,nrow(data))
cutoff2 <- rep(0,nrow(data))
cutoff2.5 <- rep(0,nrow(data))
cutoff3 <- rep(0,nrow(data))
cutoff5 <- rep(0,nrow(data))
cutoff10 <- rep(0,nrow(data))
df_curve <- data.frame(cutoff0, cutoff0.1,cutoff0.5,cutoff1, cutoff1.5, cutoff2, cutoff2.5, cutoff3,cutoff5, cutoff10)
df_curve$otuid <- rownames(data)
df_curve <- df_curve[, c(ncol(df_curve),1:(ncol(df_curve)-1))]

#_____________________________run from here________________________________________

cuts <- c(0, 0, 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1)

for (c in 2:ncol(df_curve)){
  
  otudata_tfilt <- data.frame(t(data))
  otudf <- otudata_tfilt
  otudata_tfilt <- as.data.frame(t(apply(otudata_tfilt, 1, function(x) x/sum(x))))
  rowSums(otudata_tfilt)
  #remove asv < 0.5% of reads of the total reads but keeping reads number
  for (m in 1:ncol(otudf)) {
    
    for (v in 1:nrow(otudf)) {
      
      if (otudata_tfilt[v,m] < cuts[c]) { #to account for cross contamination and barcode switching during sequencing
        
        otudf[v,m] <- 0
      }}}
  rowSums(otudf)
  colSums(otudf)
  otudf <- otudf[,colSums(otudf) > 0]
  #join metadata with asv data
  otudf <- data.frame(otudf)
  otudf$baseclear_id <- rownames(otudf)
  df10 <- merge(metadata[, c(1,2,7,8,10)], otudf, by = "baseclear_id")
  #add island number only
  splitnumber <- str_split_fixed(df10$field.nmbr, "-", 2)
  df10$island <-  splitnumber[,1]
  df10<- df10[,c(1:5,ncol(df10),6:(ncol(df10)-1))]
  
  df <- df10
  df <- melt(df, id.vars = c("baseclear_id","mastermix","registr.code","field.nmbr","taxon","island"), 
             value.name = "reads", variable.name = "Name")
  df <- df[,c(1,2,3,4,5,6,7,8)]
  df <- df[df$reads > 0, ]
  
  cutoffasv <- levels(droplevels(df$Name))
  readsum <- df %>% group_by(Name) %>% summarise(totread = sum(reads))
  
  for (x in 1:length(cutoffasv)){
    
    for(i in 1:nrow(df_curve)) {
      
      if (cutoffasv[x] == df_curve$otuid[i]) {
        
        df_curve[i,c] <- readsum$totread[x]
        
      }
    }
  }
  
}

df55 <- df_curve

rownames(df55) <- df55$otuid
df55 <- t(df55[,-1])


#plot rarefaction curve

S <- specnumber(df55) # observed number of species
raremax <- min(rowSums(df55))
Srare <- rarefy(df55, raremax)

data1 <- data
for (j in 2:ncol(data1)){
  for(i in 1:nrow(data1)) {
    if (data1[i,j] > 1) {
      
      data1[i,j] <- 1
    }
  }
  
}

colSums(data1[2:ncol(data1)])
max(colSums(data[2:ncol(data)]))

par(mfrow = c(1,2))
plot(colSums(data1[2:ncol(data1)]), colSums(data[2:ncol(data)])/1000, xlab = "Number of ASVs per sample",
     ylab = "Number of reads per sample (in thousands)", main = "Before any threshold applied")
rarecurve(df55, step = 20, col = "blue", cex = 0.6, 
          xlab = "Total number of reads", ylab = "Number of ASVs", xlim = c(0, 50000))






#--------Prep COI data - threshold 1.5---------------------------------------------------------------
#____________________1. identify samples with low ratio passed/raw reads______________________________________________________
#

df_pro <- read.csv("~/Downloads/Dataset/ASVtable-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/COI-raw-reads-table.csv", sep = ";")

df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[df_pro_t$processed_read_sum > 1000,c(135,134)]
df_raw1 <- df_raw[c(1,2)]

df_all <- merge(df_raw1, df_pro1, by = "Sample_name")
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100

summary(df_all$ratio)

#retain samples with ratio > 10% 
sample_to_analyze <- levels(droplevels(df_all$Sample_name[df_all$ratio > 10]))

#______________________________2. create working dataset with read numbers (abundance) for each sample merged with metadata________________
#load data
unoise4 <- read.csv("~/Downloads/Dataset/ASVtable-new.csv", sep = ";", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/MetadataCOIsamples20210726.csv", sep = ";")

otudata <- unoise4
#remove samples with read ratio < 10%
keep <- sample_to_analyze
otudata_tfilt<- data.frame(t(otudata[keep]))

otudf <- otudata_tfilt

otudata_tfilt <- as.data.frame(t(apply(otudata_tfilt, 1, function(x) x/sum(x))))
rowSums(otudata_tfilt)

#remove asv < 1.5% of reads of the total reads but keeping reads number
for (i in 1:ncol(otudf)) {
  
  for (x in 1:nrow(otudf)) {
    
    if (otudata_tfilt[x,i] < 0.015) { #to account for cross contamination and barcode switching during sequencing
      
      otudf[x,i] <- 0
    }
  }
}

rowSums(otudf)
colSums(otudf)

otudf <- otudf[,colSums(otudf) > 0] #from 133 asv to 51 asv

#join metadata with asv data
otudf <- data.frame(otudf)
otudf$baseclear_id <- rownames(otudf)
df_unoise4 <- merge(metadata[, c(1,2,7,8,10)], otudf, by = "baseclear_id")

#add island number only
splitnumber <- str_split_fixed(df_unoise4$field.nmbr, "-", 2)
df_unoise4$island <-  splitnumber[,1]
df_unoise4 <- df_unoise4[,c(1:5,ncol(df_unoise4),6:(ncol(df_unoise4)-1))]

#add sequences - with a melt dataset
df<- df_unoise4
Seq <- read.csv("~/Downloads/Dataset/COIsequences-new.csv", sep = ";")

#melt ASVs
df_m <- melt(df, id.vars = c("baseclear_id","mastermix","registr.code","field.nmbr","taxon","island"), 
             value.name = "reads", variable.name = "Name")
#attached sequences to ASVs
df_mseq <- merge(df_m, Seq, by = "Name")
df_mseq <- df_mseq[df_mseq$reads > 0, ]
df_mseq <- df_mseq[,c(2,3,4,5,5,6,7,1,9,8)]
write.csv(df_mseq, "~/Downloads/Dataset/COIworkingdataset-cutoff1.5.csv", row.names=FALSE)


#__________________3. create sequence text files for each species_______________________________

df <- read.csv("~/Downloads/Dataset/COIworkingdataset-cutoff1.5.csv")

#create genus and family column
df$genus <- str_split_fixed(df$taxon, " ", 2)[,1]
#keep OTUs that are in at least two specimens of a species

df2 <- df %>% group_by(genus, Name) %>% filter(n()>1) #ASV shared between at least two specimens of a species
df3 <- df %>% group_by(Name) %>% filter(n()==1) #ASV not shared between specimens
df4 <- df %>% group_by(taxon) %>% filter(n()<4) #species with 2 or less specimens
df5 <- rbind(df2, df3, df4) #combine both dataset and it should reflect the intra and inter species variation
df5 <- distinct(df5)

df5$shortcode <- str_split_fixed(df5$registr.code, ".LAB.", 2)[,2]

levels(droplevels(df5$baseclear_id))
levels(droplevels(df5$taxon))
levels(droplevels(df5$Name))
levels(droplevels(df5$registr.code))
length(levels(droplevels(df5$baseclear_id)))
length(levels(droplevels(df5$taxon)))
length(levels(droplevels(df5$Name)))
length(levels(droplevels(df5$registr.code)))

#subset species
df_sub <- df5
#new column - merge otu and species together
df_sub$asv <- paste(">",df_sub$Name,"_",df_sub$shortcode,"_",df_sub$taxon, sep = "")
df_sub$ID <- paste(df_sub$taxon,"_",df_sub$baseclear_id, sep = "")
#keep only baseclear id, asv and sequence
fortxt <- df_sub[,c("ID", "asv", "Sequence")]
fortxt <- distinct(fortxt) #remove duplicate rows

#Create loop for subsetting and write text file
fortxt <- as.data.frame(fortxt)
IDs <- unique(fortxt$ID)

for (i in 1:length(IDs)) {
  Subsetted <- subset(fortxt, ID == IDs[i])
  Subsetted <- Subsetted[,c(2,3)]
  File <- paste("~/Downloads/Dataset/cutoff_1.5/", IDs[i], ".txt", sep = "")
  write.table(Subsetted, File, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE )
}

#__________________4. create stats filtering criteria_______________________________

df_pro <- read.csv("~/Downloads/Dataset/ASVtable-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/COI-raw-reads-table.csv", sep = ";")
df <- read.csv("~/Downloads/Dataset/COIworkingdataset-cutoff1.5.csv")

#get sum of reads/sample from processed data - ASV table without cutoff
df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[,c(135,134)]
#get sum of reads/sample from raw data
df_raw1 <- df_raw[1:177,c(1,2)]
#get sum of reads/sample from ASV table with cutoff
colnames(df)[1] <- "Sample_name"
df_cut <- df %>% group_by(Sample_name) %>% summarise(readscutoff = sum(reads))

df_all <- merge(df_raw1, df_pro1, by = "Sample_name", all = TRUE)
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100
df_all <- merge(df_all, df_cut, by = "Sample_name", all = TRUE)

write.csv(df_all, "~/Downloads/Dataset/COIreadpersamplestats.csv", row.names=FALSE)


#--------prep 18S data - threshold 1.5----------------------------





#____________________1. identify samples with low ratio passed/raw reads______________________________________________________

df_pro <- read.csv("~/Downloads/Dataset/ASVtable-18S-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/18S-raw-reads-table.csv", sep = ";")

df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[df_pro_t$processed_read_sum > 1000,c(1385,1384)]
df_raw1 <- df_raw[c(1,2)]

df_all <- merge(df_raw1, df_pro1, by = "Sample_name")
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100

summary(df_all$ratio)

#retain samples with ratio > 10% 
sample_to_analyze <- levels(droplevels(df_all$Sample_name[df_all$ratio > 10]))

#______________________________2. create working dataset with read numbers (abundance) for each sample merged with metadata________________

#load data
unoise4 <- read.csv("~/Downloads/Dataset/ASVtable-18S-new.csv", sep = ";", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/Metadata18Ssamples20210726.csv", sep = ";")

otudata <- unoise4
#remove samples with read ratio < 10%
keep <- sample_to_analyze
otudata_tfilt<- data.frame(t(otudata[keep]))

otudf <- otudata_tfilt

otudata_tfilt <- as.data.frame(t(apply(otudata_tfilt, 1, function(x) x/sum(x))))
rowSums(otudata_tfilt)

#remove asv < 1.5% of reads of the total reads but keeping reads number
for (i in 1:ncol(otudf)) {
  
  for (x in 1:nrow(otudf)) {
    
    if (otudata_tfilt[x,i] < 0.015) { #to account for cross contamination and barcode switching during sequencing
      
      otudf[x,i] <- 0
    }
  }
}

rowSums(otudf)
colSums(otudf)

otudf <- otudf[,colSums(otudf) > 0] #from 1383 asv to 976 asv

#join metadata with asv data
otudf <- data.frame(otudf)
otudf$baseclear_id <- rownames(otudf)
df_unoise4 <- merge(metadata[, c(1,2,7,8,10)], otudf, by = "baseclear_id")

#add island number only
splitnumber <- str_split_fixed(df_unoise4$field.nmbr, "-", 2)
df_unoise4$island <-  splitnumber[,1]
df_unoise4 <- df_unoise4[,c(1:5,ncol(df_unoise4),6:(ncol(df_unoise4)-1))]

#add sequences - with a melt dataset
df <- df_unoise4
Seq <- read.csv("~/Downloads/Dataset/18S-sequences-new.csv", sep = ";")

#melt ASVs
df_m <- melt(df, id.vars = c("baseclear_id","mastermix","registr.code","field.nmbr","taxon","island"), 
             value.name = "reads", variable.name = "Name")
#attached sequences to ASVs
df_mseq <- merge(df_m, Seq, by = "Name")
df_mseq <- df_mseq[df_mseq$reads > 0, ]
df_mseq <- df_mseq[,c(2,3,4,5,5,6,7,1,9,8)]
write.csv(df_mseq, "~/Downloads/Dataset/18Sworkingdataset-cutoff1.5.csv", row.names=FALSE)


#__________________3. create sequence text files for each species_______________________________

df <- read.csv("~/Downloads/Dataset/18Sworkingdataset-cutoff1.5.csv")

#create genus and family column
df$genus <- str_split_fixed(df$taxon, " ", 2)[,1]
#keep OTUs that are in at least two specimens of a species

df2 <- df %>% group_by(genus, Name) %>% filter(n()>1) #ASV shared between at least two specimens of a species
df3 <- df %>% group_by(Name) %>% filter(n()==1) #ASV not shared between specimens
df4 <- df %>% group_by(taxon) %>% filter(n()<4) #species with 2 or less specimens
df5 <- rbind(df2, df3, df4) #combine both dataset and it should reflect the intra and inter species variation
df5 <- distinct(df5)

df5$shortcode <- str_split_fixed(df5$registr.code, ".LAB.", 2)[,2]

levels(droplevels(df5$baseclear_id))
levels(droplevels(df5$taxon))
levels(droplevels(df5$Name))
levels(droplevels(df5$registr.code))
length(levels(droplevels(df5$baseclear_id)))
length(levels(droplevels(df5$taxon)))
length(levels(droplevels(df5$Name)))
length(levels(droplevels(df5$registr.code)))

#subset species
df_sub <- df5
#new column - merge otu and species together
df_sub$asv <- paste(">",df_sub$Name,"_",df_sub$shortcode,"_",df_sub$taxon, sep = "")
df_sub$ID <- paste(df_sub$taxon,"_",df_sub$baseclear_id, sep = "")
#keep only baseclear id, asv and sequence
fortxt <- df_sub[,c("ID", "asv", "Sequence")]
fortxt <- distinct(fortxt) #remove duplicate rows

#Create loop for subsetting and write text file
fortxt <- as.data.frame(fortxt)
IDs <- unique(fortxt$ID)

for (i in 1:length(IDs)) {
  Subsetted <- subset(fortxt, ID == IDs[i])
  Subsetted <- Subsetted[,c(2,3)]
  File <- paste("~/Downloads/Dataset/cutoff_1.5/", IDs[i], ".txt", sep = "")
  write.table(Subsetted, File, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE )
}

#__________________4. create stats filtering criteria_______________________________

df_pro <- read.csv("~/Downloads/Dataset/ASVtable-18S-new.csv", sep = ";", row.names = 1)
df_raw <- read.csv("~/Downloads/Dataset/18S-raw-reads-table.csv", sep = ";")
df <- read.csv("~/Downloads/Dataset/18Sworkingdataset-cutoff1.5.csv")

#get sum of reads/sample from processed data - ASV table without cutoff
df_pro_t <- data.frame(t(df_pro))
df_pro_t$processed_read_sum <- rowSums(df_pro_t)
df_pro_t$Sample_name <- rownames(df_pro_t)
df_pro1 <- df_pro_t[,c(1385,1384)]
#get sum of reads/sample from raw data
df_raw1 <- df_raw[1:120,c(1,2)]
#get sum of reads/sample from ASV table with cutoff
colnames(df)[1] <- "Sample_name"
df_cut <- df %>% group_by(Sample_name) %>% summarise(readscutoff = sum(reads))

df_all <- merge(df_raw1, df_pro1, by = "Sample_name", all = TRUE)
df_all$ratio <- df_all$processed_read_sum/df_all$Number_of_read.pairs*100
df_all <- merge(df_all, df_cut, by = "Sample_name", all = TRUE)

write.csv(df_all, "~/Downloads/Dataset/18Sreadpersamplestats.csv", row.names=FALSE)


#--------18S-COI - violin chart per species for 1.5 cutoff and statistical tests-------------------------------------------------------------
df1518S <- read.csv("~/Downloads/Dataset/18Sgeneticdistancematrix1.5.csv")
df15COI <- read.csv("~/Downloads/Dataset/COIgeneticdistancematrix1.5.csv")

df21518S <- melt(df1518S, id.vars = "X", 
                 value.name = "patristic_distance", variable.name = "X2")
df21518S <- df21518S[complete.cases(df21518S),] #remove NAs
colnames(df21518S) <- c("otu_1", "otu_2", "patristic_distances") #rename columns
df21518S$cutoff <- "cutoff1.5"
df21518S$marker <- "18S"

df215COI <- melt(df15COI, id.vars = "X", 
                 value.name = "patristic_distance", variable.name = "X2")
df215COI <- df215COI[complete.cases(df215COI),] #remove NAs
colnames(df215COI) <- c("otu_1", "otu_2", "patristic_distances") #rename columns
df215COI$cutoff <- "cutoff1.5"
df215COI$marker <- "COI"

df2 <- rbind(df21518S, df215COI)

#split columns and make the matrix from geneious usable
df2$number1 <- str_split_fixed(df2$otu_1, "_", 3)[,2]
df2$number2 <- str_split_fixed(df2$otu_2, "_", 3)[,2]
df2$species1 <- str_split_fixed(df2$otu_1, "_", 3)[,3]
df2$species2 <- str_split_fixed(df2$otu_2, "_", 3)[,3]
df2$species2 <- sub("\\.", " ", df2$species2)
df2$species2 <- sub("\\.", " ", df2$species2)
df2$species1 <- sub("_", " ", df2$species1)
df2$species2 <- sub("_", " ", df2$species2)
df2$otu_2 <- sub("\\.", " ", df2$otu_2)
df2$otu_2 <- sub("\\.", " ", df2$otu_2)
df2$otu_2 <- sub("\\.", " ", df2$otu_2)
df2$species1 <- sub(" 1", "", df2$species1)
df2$species1 <- sub(" 2", "", df2$species1)
df2$species1 <- sub(" 3", "", df2$species1)
df2$species1 <- sub(" 4", "", df2$species1)
df2$species2 <- sub(" 1", "", df2$species2)
df2$species2 <- sub(" 2", "", df2$species2)
df2$species2 <- sub(" 3", "", df2$species2)
df2$species2 <- sub(" 4", "", df2$species2)
df2$otu_1 <- sub(" 1", "", df2$otu_1)
df2$otu_1 <- sub(" 2", "", df2$otu_1)
df2$otu_1 <- sub(" 3", "", df2$otu_1)
df2$otu_1 <- sub(" 4", "", df2$otu_1)
df2$otu_2 <- sub(" 1", "", df2$otu_2)
df2$otu_2 <- sub(" 2", "", df2$otu_2)
df2$otu_2 <- sub(" 3", "", df2$otu_2)
df2$otu_2 <- sub(" 4", "", df2$otu_2)
df2$asv1 <- str_split_fixed(df2$otu_1, "_", 3)[,1]
df2$asv2 <- str_split_fixed(df2$otu_2, "_", 3)[,1]
df2$asv1 <- sub(">", "", df2$asv1)
df2$asv2 <- sub("X ", "", df2$asv2)


#we want to compare different specimens, and different ASVs in the same specimen, thus remove duplicate lines
df2 <- df2[!(df2$number1 == df2$number2 & df2$asv1 == df2$asv2 & df2$marker == "COI"),]
df2 <- df2[!(df2$number1 == df2$number2 & df2$asv1 == df2$asv2 & df2$marker == "18S"),]
df2 <- distinct(df2)

df2$type <- NA
df2$type[df2$number1 == df2$number2] <- "intra"
df2$type[((df2$species1 == df2$species2) & !(df2$number1 == df2$number2))] <- "inter"
df2$type[!(df2$species1 == df2$species2)] <- "between"

species_for_comparison <- c("Nummulites venosus","Amphisorus SpL","Amphisorus SpS","Parasorites sp","Marginopora vertebralis",
                            "Operculina ammonoides","Amphistegina radiata")

df3 <- df2[df2$species1 %in% species_for_comparison &
             df2$species2 %in% species_for_comparison &
             df2$type %in% c("intra", "inter"),]

#function to produce summary statistic - mean +/- standard deviation
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#per species
ggplot(df3, aes(marker, patristic_distances)) + 
  geom_violin(aes(color = marker, fill = marker), alpha = 0.25) + 
  stat_summary(fun.data=data_summary) +
  scale_fill_manual(values = c("18S" = "#CC0033", "COI" = "#0066CC")) +
  scale_color_manual(values = c("18S" = "#CC0033", "COI" = "#0066CC")) +
  theme_bw() + ylim(0,1) + labs(y = "Patristic distance") + facet_grid(type~species1) +
  #coord_flip() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x = element_blank(), legend.position = "bottom")

#statistical test 

Species <- species_for_comparison
Intraleve <- rep(0, length(species_for_comparison))
Intraflig <- rep(0, length(species_for_comparison))
IntraleveF <- rep(0, length(species_for_comparison))
IntrafligC <- rep(0, length(species_for_comparison))
Interleve <- rep(0, length(species_for_comparison))
Interflig <- rep(0, length(species_for_comparison))
InterleveF <- rep(0, length(species_for_comparison))
InterfligC <- rep(0, length(species_for_comparison))
testres <- data.frame(Species, IntraleveF, Intraleve, IntrafligC, Intraflig, InterleveF, Interleve, InterfligC, Interflig)

for(i in 1:length(species_for_comparison)) {
  
  df4 <- df3[df3$species1 %in% species_for_comparison[i] &
               df3$type %in% "intra",]
  df4$marker <- as.factor(df4$marker)
  #two options: Levene's test (better than ftest) and fligner-killeen test (better for non-normally distribution)
  
  if (length(unique(df4$marker)) < 2) {
    
    testres$Intraleve[i] <- NA
    testres$Intraflig[i] <- NA
    
    df4 <- df3[df3$species1 %in% species_for_comparison[i] &
                 df3$type %in% "inter",]
    df4$marker <- as.factor(df4$marker)
    
    if (length(unique(df4$marker)) < 2) {
      
      testres$Interleve[i] <- NA
      testres$Interflig[i] <- NA
      
      break
    }
    
    #two options: Levene's test (better than ftest) and fligner-killeen test (better for non-normally distribution)
    leve <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,3]
    leveF <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,2]
    flig <- fligner.test(patristic_distances ~ marker, data = df4)[3]
    fligC <- fligner.test(patristic_distances ~ marker, data = df4)[1]
    
    testres$Interleve[i] <- leve
    testres$Interflig[i] <- flig
    testres$InterleveF[i] <- leveF
    testres$InterfligC[i] <- fligC
    
  } else {
    
    leve <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,3]
    leveF <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,2]
    flig <- fligner.test(patristic_distances ~ marker, data = df4)[3]
    fligC <- fligner.test(patristic_distances ~ marker, data = df4)[1]
    
    testres$Intraleve[i] <- leve
    testres$Intraflig[i] <- flig
    testres$IntraleveF[i] <- leveF
    testres$IntrafligC[i] <- fligC
    
    df4 <- df3[df3$species1 %in% species_for_comparison[i] &
                 df3$type %in% "inter",]
    df4$marker <- as.factor(df4$marker)
    #two options: Levene's test (better than ftest) and fligner-killeen test (better for non-normally distribution)
    leve <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,3]
    leveF <- leveneTest(y = df4$patristic_distances, group = df4$marker, center=mean)[1,2]
    flig <- fligner.test(patristic_distances ~ marker, data = df4)[3]
    fligC <- fligner.test(patristic_distances ~ marker, data = df4)[1]
    
    testres$Interleve[i] <- leve
    testres$Interflig[i] <- flig
    testres$InterleveF[i] <- leveF
    testres$InterfligC[i] <- fligC
  }
  
  
}


testres <- apply(testres,2,as.character)

write.csv(testres, "~/Downloads/Dataset/statisticaltestsres.csv", row.names = FALSE)



#---------18S-COI - barchart of number of ASVs per samples --------------------------------------------
dfCOI <- read.csv("~/Downloads/Dataset/COIworkingdataset-cutoff1.5.csv")
#create genus and family column
dfCOI$genus <- str_split_fixed(dfCOI$taxon, " ", 2)[,1]
#keep OTUs that are in at least two specimens of a species
dfCOI2 <- dfCOI %>% group_by(genus, Name) %>% filter(n()>1) #ASV shared between at least two specimens of a species
dfCOI3 <- dfCOI %>% group_by(Name) %>% filter(n()==1) #ASV not shared between specimens
dfCOI4 <- dfCOI %>% group_by(taxon) %>% filter(n()<4) #species with 2 or less specimens
dfCOI5 <- rbind(dfCOI2, dfCOI3, dfCOI4) #combine both dataset and it should reflect the intra and inter species variation
dfCOI5 <- distinct(dfCOI5)
dfCOI5$shortcode <- str_split_fixed(dfCOI5$registr.code, ".LAB.", 2)[,2]
dfCOI5$marker <- "COI"
levels(droplevels(dfCOI5$baseclear_id))
levels(droplevels(dfCOI5$taxon))
levels(droplevels(dfCOI5$Name))
levels(droplevels(dfCOI5$registr.code))

df18S <- read.csv("~/Downloads/Dataset/18Sworkingdataset-cutoff1.5.csv")
#create genus and family column
df18S$genus <- str_split_fixed(df18S$taxon, " ", 2)[,1]
#keep OTUs that are in at least two specimens of a species
df18S2 <- df18S %>% group_by(genus, Name) %>% filter(n()>1) #ASV shared between at least two specimens of a species
df18S3 <- df18S %>% group_by(Name) %>% filter(n()==1) #ASV not shared between specimens
df18S4 <- df18S %>% group_by(taxon) %>% filter(n()<4) #species with 2 or less specimens
df18S5 <- rbind(df18S2, df18S3, df18S4) #combine both dataset and it should reflect the intra and inter species variation
df18S5 <- distinct(df18S5)
df18S5$shortcode <- str_split_fixed(df18S5$registr.code, ".LAB.", 2)[,2]
df18S5$marker <- "18S"
levels(droplevels(df18S5$baseclear_id))
levels(droplevels(df18S5$taxon))
levels(droplevels(df18S5$Name))
levels(droplevels(df18S5$registr.code))

df <- rbind(df18S5, dfCOI5)

species_for_comparison <- c("Nummulites venosus","Amphisorus SpL","Amphisorus SpS","Parasorites sp","Marginopora vertebralis",
                            "Operculina ammonoides","Amphistegina radiata")

df3 <- df[df$taxon %in% species_for_comparison,]==

df3$asvcount <- 1

df3 <- with(df3, df3[order(taxon,registr.code),])
order <- unique(df3$registr.code)
df3$registr.code <- factor(df3$registr.code, levels = order)

df4 <- df3[,c("registr.code", "asvcount", "marker", "Name", "taxon")]
df4 <- distinct(df4)

df5 <- df4 %>% group_by(taxon, registr.code, marker) %>% summarise(totalasv = sum(asvcount)) #compile how many asvs per sample

ggplot(df5, aes(x = marker, y = totalasv)) + 
  geom_boxplot(aes(color = marker, fill = marker), alpha = 0.25) + 
  theme_bw() + 
  scale_fill_manual(values = c("18S" = "#CC0033", "COI" = "#0066CC")) +
  scale_color_manual(values = c("18S" = "#CC0033", "COI" = "#0066CC")) +
  facet_grid(~taxon) +
  ggtitle("Number of ASVs per species") +
  xlab("Species") + ylab("Number of ASVs")






#--------18S-COI - Phylogenetic tree plot---------------------------------------------


#______________________COI tree arrangement__________________________________________________________________________________
tree <- read.nexus("~/Downloads/Dataset/COItree-rerooted-boot.txt")
tree
tree <- midpoint.root(tree)

tree <- groupClade(tree, c(346, 398, 355, 356, 482, 264, 255, 257, 258, 259))
nodes <- c(346, 398, 355, 356, 482, 264, 255, 257, 258, 259)
boot <- tree$node.label
tips <- rep("", 252)
bootstrap <- c(tips, boot)

#to make node colors below 75 bootstrap support grey.
Col <- rep("black", 251)
booty <- as.data.frame(cbind(boot, Col))
booty$boot <- as.numeric(as.character(booty$boot))
booty$Col <- as.character(booty$Col)
booty$boot[is.na(booty$boot)] <- 100
booty$Col[booty$boot < 75] <- "red"

species <- as.data.frame(tree[["tip.label"]])
species$sp <- str_split_fixed(species$`tree[["tip.label"]]`, "_", 4)[,4]
species$sp <- sub("_1", "", species$sp)
species$sp <- sub("_2", "", species$sp)
species$sp <- sub("_3", "", species$sp)
species$sp <- sub("_4", "", species$sp)
species$fam[species$sp %in% c("Sorites_sp","Marginopora_vertebralis","Parasorites_sp", "Amphisorus_SpS","Amphisorus_SpL")] <- "Soritidae"
species$fam[species$sp %in% c("Calcarina_sp2","Calcarina_sp1", "Neorotalia_calcar", "Neorotalia_gaimardi")] <- "Calcarinidae"
species$fam[species$sp %in% c("Amphistegina_papillosa", "Amphistegina_lessonii", "Amphistegina_radiata")] <- "Amphisteginidae"
species$fam[species$sp %in% c("Peneroplis_planatus")] <- "Peneroplidae"
species$fam[species$sp %in% c("Alveolinella_quoyi", "Borelis_schlumbergeri")] <- "Alveolinidae"
species$fam[species$sp %in% c("Nummulites_venosus", "Operculina_sp1", "Operculinella_cummingi","Operculina_ammonoides", "Operculina_LKI27type", "Heterostegina_depressa", "Operculina_complanata")] <- "Nummulitidae"
row.names(species) <- NULL
species$shap[species$sp %in% c("Sorites_sp","Peneroplis_planatus", "Alveolinella_quoyi", "Nummulites_venosus","Amphistegina_radiata","Calcarina_sp2")] <- "A"
species$shap[species$sp %in% c("Parasorites_sp", "Borelis_schlumbergeri", "Operculina_ammonoides","Calcarina_sp1", "Operculina_complanata","Amphistegina_papillosa")] <- "B"
species$shap[species$sp %in% c("Marginopora_vertebralis","Neorotalia_gaimardi","Operculina_sp1","Amphistegina_lessonii")] <- "C"
species$shap[species$sp %in% c("Amphisorus_SpS","Neorotalia_calcar","Operculina_LKI27type")] <- "D"
species$shap[species$sp %in% c("Amphisorus_SpL","Operculinella_cummingi", "Operculina_complanata")] <- "E"
species$shap[species$sp %in% c("Heterostegina_depressa")] <- "F"
species$shap[species$sp %in% c("Operculinella_cummingi")] <- "G"

p <- ggtree(tree, layout = "circular", size = 1) + 
  geom_tiplab(align = TRUE, linesize=0.5) 
#  + geom_text(aes(label=bootstrap), color = "red") 
  
p <- p + geom_nodepoint(color=booty$Col, size=3)

p <- p %<+% species + geom_tippoint(aes(color=fam, shape=shap, x=0.68), size=4) +
  scale_color_manual(values=c("Alveolinidae" = "#FF9933", "Amphisteginidae"="#CC33CC", "Calcarinidae"="#CC0033", "Nummulitidae"="#0066CC","Peneroplidae"="#00CC33","Soritidae"="#FFCC00")) +
  scale_shape_manual(values=c(15,16,17,0,1,2,8)) 

plot(p)




#______________________18S tree arrangement__________________________________________________________________________________
tree2 <- read.nexus("~/Downloads/Dataset/18Stree-rerooted-boot.txt")
tree2
tree2 <- midpoint.root(tree2)

boot2 <- tree2$node.label
tips2 <- rep("", 504)
bootstrap2 <- c(tips2, boot2)

#to make node colors below 75 bootstrap support grey.
Col2 <- rep("black", 503)
booty2 <- as.data.frame(cbind(boot2, Col2))
booty2$boot2 <- as.numeric(as.character(booty2$boot2))
booty2$Col2 <- as.character(booty2$Col2)
booty2$boot2[is.na(booty2$boot2)] <- 100
booty2$Col2[booty2$boot2 < 75] <- "red"

species <- as.data.frame(tree2[["tip.label"]])
species$sp <- str_split_fixed(species$`tree2[["tip.label"]]`, "_", 4)[,4]
species$sp <- sub("_1", "", species$sp)
species$sp <- sub("_2", "", species$sp)
species$sp <- sub("_3", "", species$sp)
species$sp <- sub("_4", "", species$sp)
species$fam[species$sp %in% c("Sorites_sp","Marginopora_vertebralis","Parasorites_sp", "Amphisorus_SpS","Amphisorus_SpL")] <- "Soritidae"
species$fam[species$sp %in% c("Calcarina_sp2","Calcarina_sp1", "Neorotalia_calcar", "Neorotalia_gaimardi")] <- "Calcarinidae"
species$fam[species$sp %in% c("Amphistegina_radiata")] <- "Amphisteginidae"
species$fam[species$sp %in% c("Alveolinella_quoyi", "Borelis_schlumbergeri")] <- "Alveolinidae"
species$fam[species$sp %in% c("Nummulites_venosus", "Operculinella_cummingi","Operculina_ammonoides", "Operculina_LKI27type","Operculina_complanata")] <- "Nummulitidae"
row.names(species) <- NULL
species$shap[species$sp %in% c("Sorites_sp","Alveolinella_quoyi", "Nummulites_venosus","Amphistegina_radiata","Calcarina_sp2")] <- "A"
species$shap[species$sp %in% c("Parasorites_sp", "Borelis_schlumbergeri", "Operculina_ammonoides","Calcarina_sp1", "Operculina_complanata")] <- "B"
species$shap[species$sp %in% c("Marginopora_vertebralis","Neorotalia_gaimardi","Operculina_sp1")] <- "C"
species$shap[species$sp %in% c("Amphisorus_SpS","Neorotalia_calcar","Operculina_LKI27type")] <- "D"
species$shap[species$sp %in% c("Amphisorus_SpL","Operculinella_cummingi", "Operculina_complanata")] <- "E"
species$shap[species$sp %in% c("Operculinella_cummingi")] <- "G"

p <- ggtree(tree2, layout = "circular", size = 1) + 
  geom_tiplab(align = TRUE, linesize=0.5) 
#  + geom_text(aes(label=bootstrap2), color = "red")

p <- p + geom_nodepoint(color=booty2$Col2, size=3)

p <- p %<+% species + geom_tippoint(aes(color=fam, shape=shap, x=1.29), size=2) +
  scale_color_manual(values=c("Alveolinidae" = "#FF9933", "Amphisteginidae"="#CC33CC", "Calcarinidae"="#CC0033", "Nummulitidae"="#0066CC","Soritidae"="#FFCC00")) +
  scale_shape_manual(values=c(15,16,17,0,1,8))

plot(p)


