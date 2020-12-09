from sklearn import svm, tree, linear_model
import numpy as np
from graphviz import Source
from IPython.display import SVG

import csv

#Script to make a support vector machine and decision tree which predicts drug resistance from genetic profiles

#Note code requires pre-processing performed in R ahead of time

#Open the data frame with both the mutations and the resistance
print("Opening data files")
with open('mutationsIntegerDataFrame.csv', newline='') as f:
    reader = csv.reader(f)
    mutationsDataFrame = list(reader)


with open('resistanceIntegerDataFrame.csv', newline='') as f2:
    reader = csv.reader(f2)
    resistanceDataFrame = list(reader)


#Convert data frames to integer form
print("Converting data to integer form")
for i in range(len(mutationsDataFrame)):
     mutationsDataFrame[i] = list(map(int, mutationsDataFrame[i]))

for i in range(len(resistanceDataFrame)):
     resistanceDataFrame[i] = int(resistanceDataFrame[i][0])


#Make the training and testing datasets by dividing the data in half
print("Constructing training and testing datasets")
mutTrainingData = mutationsDataFrame[1:int(len(mutationsDataFrame)/2)]
mutTestingData = mutationsDataFrame[int(len(mutationsDataFrame)/2)+1 :
                                    len(mutationsDataFrame)]

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
                    clf.fit(mutTrainingData,resisTrainingData)

                    trainDatPredictions = list(clf.predict(mutTrainingData))
                    testDatPredictions = list(clf.predict(mutTestingData))

                    correctList = find_accuracy(trainDatPredictions, testDatPredictions)
                    datList.append(correctList)

                    correctList.insert(0,c)
                    correctList.insert(0,deg)
                              
                    csvwriter.writerow(correctList)




#Train using tree analysis
print("\n\nTraining the decision tree\n")
clf2 = tree.DecisionTreeClassifier(min_impurity_decrease = 0.01)
clf2 = clf2.fit(mutTrainingData,resisTrainingData)



trainDatPredictions = list(clf2.predict(mutTrainingData))
testDatPredictions = list(clf2.predict(mutTestingData))

print(find_accuracy(trainDatPredictions, testDatPredictions))

feat_names = ["embB","gyrA","inhA","katG","PncA","rpoB","rpsL","rrs"]
class_names = ["MDR non XDR", "Mono DR", "Poly DR", "Sensitive", "XDR"]
graph = Source(tree.export_graphviz(clf2, out_file = None, feature_names = feat_names, class_names = class_names, filled = True, impurity = False))
graph.format = 'png'
graph.render('dtree_render')



