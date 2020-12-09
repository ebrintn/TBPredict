#Program to clean genomics data for Tuberculosis preliminary investigation. Note the code requires 
#access to the TB portals genomic data which cannot be shared but can be accessed with permission
#from TB portals


#First load in genomics data and patient data into the program 
library(readr)
TB_Portals_Genomics_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals_Genomics_20201021.csv")
TB_Portals_Patient_Cases_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Patient Cases_20201021.csv")


#Remove dubplicated data
TB_Portals_Genomics_20201021_unique<-TB_Portals_Genomics_20201021[!duplicated(TB_Portals_Genomics_20201021$condition_id), ]
TB_Portals_Genomics_20201021_unique <- TB_Portals_Genomics_20201021_unique[order(TB_Portals_Genomics_20201021_unique$condition_id),]

ids <- TB_Portals_Genomics_20201021_unique$condition_id


#Find the appropriate patient samples for the samples with genomic information
neededSamples <- TB_Portals_Patient_Cases_20201021$condition_id%in%ids
availableSamples <- TB_Portals_Patient_Cases_20201021$condition_id[neededSamples]
sampleResistanceType <- TB_Portals_Patient_Cases_20201021$type_of_resistance[neededSamples]

sampleMatrix <- matrix(cbind(availableSamples, sampleResistanceType),ncol = 2)
sampleMatrix <- sampleMatrix[order(sampleMatrix[,1]),]
sampleResistanceType <- sampleMatrix[,2]

#Record genomic resistance type as integers
integerSampleResistanceType <- c()
for(sample in sampleResistanceType){
  if(sample == "MDR non XDR"){
    resisType = 0
  } else if (sample == "Mono DR"){
    resisType = 1
  } else if (sample == "Poly DR"){
    resisType = 2
  } else if (sample == "Sensitive"){
    resisType = 3
  } else if (sample == "XDR"){
    resisType = 4
  }
  integerSampleResistanceType = c(integerSampleResistanceType, resisType)
}




#Find the required SNPs, record as integers
prevalentSNPS <- TB_Portals_Genomics_20201021_unique$gene_snp_mutations

embBSNPS <- as.integer(as.logical(grepl("embB", prevalentSNPS)))
gyrASNPS <- as.integer(as.logical(grepl("gyrA", prevalentSNPS)))
inhAProSNPS <- as.integer(as.logical(grepl("inhA-Pro", prevalentSNPS)))
katGSNPS <- as.integer(as.logical(grepl("katG", prevalentSNPS)))
pncASNPS <- as.integer(as.logical(grepl("pncA", prevalentSNPS)))
rpoBSNPS <- as.integer(as.logical(grepl("rpoB", prevalentSNPS)))
rpsLSNPS <- as.integer(as.logical(grepl("rpsL", prevalentSNPS)))
rrsSNPS <- as.integer(as.logical(grepl("rrs", prevalentSNPS)))


resistanceMatrix <- data.frame(cbind( embBSNPS, gyrASNPS, inhAProSNPS, katGSNPS, pncASNPS, rpoBSNPS, rpsLSNPS, rrsSNPS))


#Look at PCA of all data samples
write.csv(resistanceMatrix,"mutationsIntegerDataFrame.csv", row.names = FALSE)
write.csv(integerSampleResistanceType, "resistanceIntegerDataFrame.csv", row.names = F)



#Convert via PCA
library(ggfortify)
resistance.pca <- prcomp(resistanceMatrix, center = TRUE, scale = TRUE)
colnames(sampleMatrix) <- c("id","phenotype")
autoplot(resistance.pca, data = sampleMatrix, colour = "phenotype")
write.csv(resistance.pca$x[,1:2], "pcaResistanceDataFrame.csv", row.names = F)

#PCA of training and testing sample
training_sample <- resistanceMatrix[1:(nrow(resistanceMatrix)/2),]
testing_sample <- resistanceMatrix[(nrow(resistanceMatrix)/2)+1:nrow(resistanceMatrix),]
resistance.pca.training <- prcomp(training_sample, center = TRUE, scale = TRUE)
resistance.pca.testing <- predict(resistance.pca.training, newdata = testing_sample)
autoplot(resistance.pca.training, data = sampleMatrix[1:(nrow(sampleMatrix)/2),], colour = "phenotype")

recombined <- rbind(resistance.pca.training$x[,1:2],resistance.pca.testing[,1:2])
write.csv(recombined, "pcaMutationsDataFrameTrainingTesting.csv", row.names = F)
