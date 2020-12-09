library(readr)

#Note some preprocessing was done in excel on the files to fill in with average data


#Set up matrices of data
TB_Portals_Patient_Cases_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Patient Cases_20201021.csv")
TB_Portals_Biochemistry_20201021_Filtered <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Biochemistry_20201021 - Filtered.csv")


biochemistry_matrix <- as.matrix(TB_Portals_Biochemistry_20201021_Filtered[!duplicated(TB_Portals_Biochemistry_20201021_Filtered$patient_id), ])
biochemistry_ordered <- biochemistry_matrix[order(biochemistry_matrix[,1]),]

ids <- biochemistry_matrix[,1]


neededSamples <- TB_Portals_Patient_Cases_20201021$patient_id%in%ids
availableSamples <- TB_Portals_Patient_Cases_20201021$patient_id[neededSamples]
sampleResistanceType <- TB_Portals_Patient_Cases_20201021$type_of_resistance[neededSamples]
sampleMatrix <- matrix(cbind(availableSamples, sampleResistanceType),ncol = 2)
sampleMatrix <- sampleMatrix[order(sampleMatrix[,1]),]
sampleMatrix <- sampleMatrix[!duplicated(sampleMatrix[,1]),] 
sampleMatrix <- sampleMatrix[,2]

#Record genomic resistance type as integers
integerSampleResistanceType <- c()
for(sample in sampleMatrix){
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
  else{print("This should not be printing")}
  integerSampleResistanceType = c(integerSampleResistanceType, resisType)
}


#Write required matrices for SVM and decision tree analysis
biochemistry <- biochemistry_matrix[,2:ncol(biochemistry_matrix)]
biochem <- matrix(sapply(biochemistry, as.numeric), nrow(biochemistry), ncol(biochemistry))
colnames(biochem)<- colnames(biochemistry)
write.csv(biochemistry,"biochemistryDataFrame.csv", row.names = FALSE)
write.csv(integerSampleResistanceType, "resistanceIntegerDataFrame.csv", row.names = F)

#PCA of both training and testing
biochem.pca <- prcomp(biochem, center = T, scale = T)
sampleMatrix <- as.matrix(sampleMatrix)
colnames(sampleMatrix)<-"phenotype"
autoplot(biochem.pca, data = sampleMatrix, colour = "phenotype")


#PCA of training and testing sample
training_sample <- biochem[1:(nrow(biochem)/2),]
testing_sample <- biochem[((nrow(biochem)/2)+1):nrow(biochem),]
biochem.pca.training <- prcomp(training_sample, center = TRUE, scale = TRUE)
biochem.pca.testing <- predict(biochem.pca.training, newdata = testing_sample)
half_sample_matrix <- as.matrix(sampleMatrix[1:(nrow(sampleMatrix)/2)])
colnames(half_sample_matrix)<-"phenotype"
autoplot(biochem.pca.training, data = half_sample_matrix, colour = "phenotype")

recombined <- rbind(biochem.pca.training$x[,1:2],biochem.pca.testing[,1:2])
write.csv(recombined, "pcaBiochemicalDataFrameTrainingTesting.csv", row.names = F)
