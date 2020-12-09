from sklearn import svm, tree, linear_model
import numpy as np
from graphviz import Source
from IPython.display import SVG

import csv

#Script to make a support vector machine and decision tree which predicts drug resistance from genetic profiles

#Note code requires pre-processing performed in R ahead of time

#Open the data frame with both the mutations/genetic data and the resistance
print("Opening data files")
with open('genomicBiochemicalDataFrame.csv', newline='') as f:
    reader = csv.reader(f)
    genomicBiochemicalDataFrame = list(reader)


with open('resistanceIntegerDataFrame.csv', newline='') as f2:
    reader = csv.reader(f2)
    resistanceDataFrame = list(reader)


#Convert data frames to integer form
print("Converting data to integer form")
for i in range(len(genomicBiochemicalDataFrame)):
     genomicBiochemicalDataFrame[i] = list(map(float, genomicBiochemicalDataFrame[i]))

for i in range(len(resistanceDataFrame)):
     resistanceDataFrame[i] = int(resistanceDataFrame[i][0])


#Make the training and testing datasets by dividing the data in half
print("Constructing training and testing datasets")
genomicBiochemicaTrainingData = genomicBiochemicalDataFrame[1:int(len(genomicBiochemicalDataFrame)/2)]
genomicBiochemicaTestingData = genomicBiochemicalDataFrame[int(len(genomicBiochemicalDataFrame)/2)+1 :
                                    len(genomicBiochemicalDataFrame)]

resisTrainingData = resistanceDataFrame[1:int(len(resistanceDataFrame)/2)]
resisTestingData = resistanceDataFrame[int(len(resistanceDataFrame)/2)+1 :
                                     len(resistanceDataFrame)]


#Function to assess the accuracy of the data
def find_accuracy(predictions_training, predictions_testing,
                  training_data = resisTrainingData ,testing_data = resisTestingData):

     numCorrect = 0

     for i in range (len(training_data)):
         if (training_data[i] == predictions_training[i]):
            numCorrect += 1

     percentCorrectTest = numCorrect/len(training_data)


     numCorrect = 0

     for i in range (len(testing_data)):
         if (testing_data[i] == predictions_testing[i]):
            numCorrect += 1

     percentCorrectTrain = numCorrect/len(testing_data)

     return ([percentCorrectTest,percentCorrectTrain])




#Train using a support vector machine
for k in ('rbf','linear','poly','sigmoid'):
     datList = []

     with open('PredictionsOfEachKernel/'+ k+'.csv', 'w', newline='') as csvfile:
          csvwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
          csvwriter.writerow(["Degree","C","Training Data Validation","Testing Data Validation"])
          for c in (1,10,100,1000):
               for deg in (3,5,7,9):
                    
                    print("Training the SVM - kernel = ",k," - C = ",c," - deg = ",deg,"\n" )
                    clf = svm.SVC(kernel = k, C =c, degree = deg)
                    clf.fit(genomicBiochemicaTrainingData,resisTrainingData)

                    trainDatPredictions = list(clf.predict(genomicBiochemicaTrainingData))
                    testDatPredictions = list(clf.predict(genomicBiochemicaTestingData))

                    correctList = find_accuracy(trainDatPredictions, testDatPredictions)
                    datList.append(correctList)

                    correctList.insert(0,c)
                    correctList.insert(0,deg)
                              
                    csvwriter.writerow(correctList)




#Train using tree analysis
print("\n\nTraining the decision tree\n")
clf2 = tree.DecisionTreeClassifier(min_impurity_decrease = 0.01)
clf2 = clf2.fit(genomicBiochemicaTrainingData,resisTrainingData)



trainDatPredictions = list(clf2.predict(genomicBiochemicaTrainingData))
testDatPredictions = list(clf2.predict(genomicBiochemicaTestingData))

print(find_accuracy(trainDatPredictions, testDatPredictions))

#Make an a figure of the tree
feat_names = [
    "Thyroid Stimulating Hormone (miu/l)","Total Protein (g/l)","Potassium (mmol/l)",
    "Aspartate Aminotransferase (u/l)","Hepatitis B Surface Antigen","Hepatitis C Virus",
    "Glucose: Diabetic Post-Meal (mmol/l)","Albumin (g/l)","Urea (mmol/l)",
    "Alkaline phosphatase (u/l)","Total Bilirubin (µmol/l)","Lipase (u/l)",
    "International Normalized Ratio","Glucose (mmol/l)","Sodium (mmol/l)",
    "Thrombin Time (seconds)","Prothrombin Time (%)","Creatinine (µmol/l)",
    "Alanine Aminotransferase (u/l)","embBSNPS","gyrASNPS","inhAProSNPS",
    "katGSNPS","pncASNPS","rpoBSNPS","rpsLSNPS","rrsSNPS"
    ]
class_names = ["MDR non XDR", "Mono DR", "Poly DR", "Sensitive", "XDR"]
graph = Source(tree.export_graphviz(clf2, out_file = None, feature_names = feat_names, class_names = class_names, filled = True, impurity = False))
graph.format = 'png'
graph.render('dtree_render')



