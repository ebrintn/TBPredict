library(readr)

#Note data is available through TB portals, some preprocessing of the Biochemistry data was performed in excel
#Load required data
TB_Portals_Genomics_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals_Genomics_20201021.csv")
TB_Portals_Biochemistry_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Biochemistry_20201021.csv")
TB_Portals_Biochemistry_20201021_Filtered <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Biochemistry_20201021 - Filtered.csv")
TB_Portals_Patient_Cases_20201021 <- read_csv("~/4th_Year/Semester_1/MDSC_523/Project/Locker/TB Portals Published data_20201021/TB Portals Published data_20201021/TB Portals Patient Cases_20201021.csv")


#Remove dubplicated data
TB_Portals_Genomics_20201021<-TB_Portals_Genomics_20201021[!duplicated(TB_Portals_Genomics_20201021$condition_id), ]
TB_Portals_Biochemistry_20201021<-TB_Portals_Biochemistry_20201021[!duplicated(TB_Portals_Biochemistry_20201021$condition_id), ]
TB_Portals_Biochemistry_20201021_Filtered<-TB_Portals_Biochemistry_20201021_Filtered[!duplicated(TB_Portals_Biochemistry_20201021_Filtered$patient_id), ]
TB_Portals_Biochemistry_20201021_Filtered$condition_id <- TB_Portals_Biochemistry_20201021$condition_id

#Get the condition ids of the samples that have both genomic and biochemical data
biochem_condition_ids <- TB_Portals_Biochemistry_20201021$condition_id[TB_Portals_Biochemistry_20201021$patient_id %in%                                                                       TB_Portals_Biochemistry_20201021_Filtered$patient_id]
genomic_condition_ids <- TB_Portals_Genomics_20201021$condition_id
ids_in_both <- biochem_condition_ids[biochem_condition_ids %in% genomic_condition_ids]



#Find required data for genomic, phenotypic and biochemistry data, order data in same order
genomic <- TB_Portals_Genomics_20201021[TB_Portals_Genomics_20201021$condition_id %in% ids_in_both,]
biochemical <- TB_Portals_Biochemistry_20201021_Filtered[TB_Portals_Biochemistry_20201021_Filtered$condition_id%in%ids_in_both,]
phenotypic <- TB_Portals_Patient_Cases_20201021[TB_Portals_Patient_Cases_20201021$condition_id%in%ids_in_both,]

genomic <- genomic[order(genomic$condition_id),]
biochemical <- biochemical[order(biochemical$condition_id),]
phenotypic <- phenotypic[order(phenotypic$condition_id),]


#Remove missing genomic data
genomic <- genomic$gene_snp_mutations
snps_to_keep <- !is.na(genomic)
genomic <- genomic[snps_to_keep]
biochemical <- biochemical[snps_to_keep, 2:20]
phenotypic <- phenotypic$type_of_resistance[snps_to_keep]

print(paste("There are " , length(genomic) , "samples in both biochem and genomic data"))

#Record phenotypic presentation as integers
integerphenotype <- c()
for(sample in phenotypic){
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
  integerphenotype = c(integerphenotype, resisType)
}

#Record genomic presentation as integers
embBSNPS <- as.integer(as.logical(grepl("embB", genomic)))
gyrASNPS <- as.integer(as.logical(grepl("gyrA", genomic)))
inhAProSNPS <- as.integer(as.logical(grepl("inhA-Pro", genomic)))
katGSNPS <- as.integer(as.logical(grepl("katG", genomic)))
pncASNPS <- as.integer(as.logical(grepl("pncA", genomic)))
rpoBSNPS <- as.integer(as.logical(grepl("rpoB", genomic)))
rpsLSNPS <- as.integer(as.logical(grepl("rpsL", genomic)))
rrsSNPS <- as.integer(as.logical(grepl("rrs", genomic)))


genomicMatrix <- data.frame(cbind( embBSNPS, gyrASNPS, inhAProSNPS, katGSNPS, pncASNPS, rpoBSNPS, rpsLSNPS, rrsSNPS))


#Make a matrix of both the genomic and biochemical data
combinedMatrix <- cbind(biochemical, genomicMatrix)


#Write out the data that has been preprocessed 
write.csv(combinedMatrix,"genomicBiochemicalDataFrame.csv", row.names = FALSE)
write.csv(integerphenotype, "resistanceIntegerDataFrame.csv", row.names = F)


#Perform PCA on the whole combined Matrix
library(ggfortify)
combined.pca <- prcomp(combinedMatrix, center = TRUE, scale = TRUE)
phenotypic <- as.matrix(phenotypic)
colnames(phenotypic) <- "phenotype"
autoplot(combined.pca, data = phenotypic, colour = "phenotype")


#Perform PCA on just the training data and fit the testing data
training_sample <- combinedMatrix[1:(nrow(combinedMatrix)/2),]
testing_sample <- combinedMatrix[(nrow(combinedMatrix)/2)+1:nrow(combinedMatrix),]
combined.pca.training <- prcomp(training_sample, center = TRUE, scale = TRUE)
combined.pca.testing <- predict(combined.pca.training, newdata = testing_sample)
phenotypic.training <- as.matrix(phenotypic[1:(nrow(combinedMatrix)/2),])
colnames(phenotypic.training) <- "phenotype"
autoplot(combined.pca.training, data = phenotypic.training, colour = "phenotype")

#Write the PCA matrix
recombined <- rbind(combined.pca.training$x[,1:2],combined.pca.testing[,1:2])
write.csv(recombined, "pcaGenomicBiochemicalDataFrameTrainingTesting.csv", row.names = F)
