from sklearn import svm, tree, linear_model, datasets
import numpy as np
import matplotlib.pyplot as plt
import csv

#Script to classify genetic data using SVMs and decision trees from first 2 dimensions of PCA
#Plotting of this decision function also occurs in this script
#Note: Plotting function is adapted from https://scikit-learn.org/stable/auto_examples/svm/plot_iris_svc.html

#Open the data frame with both the mutations and the resistance
print("Opening data files")
with open('pcaResistanceDataFrameTrainingTesting.csv', newline='') as f:
    reader = csv.reader(f)
    mutationsDataFrame = list(reader)


with open('resistanceIntegerDataFrame2.csv', newline='') as f2:
    reader = csv.reader(f2)
    resistanceDataFrame = list(reader)

#Convert data frames to integer form
print("Converting data to integer form")
for i in range(len(mutationsDataFrame)):
     mutationsDataFrame[i] = list(map(float, mutationsDataFrame[i]))

for i in range(len(resistanceDataFrame)):
     resistanceDataFrame[i] = int(resistanceDataFrame[i][0])

print(mutationsDataFrame)
print(resistanceDataFrame)

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
     #print("Percent Correct Training: ", percentCorrect)


     numCorrect = 0

     for i in range (len(testing_data)):
         if (testing_data[i] == predictions_testing[i]):
            numCorrect += 1

     percentCorrectTrain = numCorrect/len(testing_data)
     #print("Percent Correct Testing: ", percentCorrect)

     return ([percentCorrectTest,percentCorrectTrain])




#Train using a support vector machine
for k in ('rbf','linear','poly','sigmoid'):
     datList = []

     with open('PredictionsOfEachKernel/'+ k+'_pca.csv', 'w', newline='') as csvfile:
          csvwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
          csvwriter.writerow(["Degree","C","Training Data Validation","Testing Data Validation"])
          for c in (1,10,100):
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
clf2 = tree.DecisionTreeClassifier()
clf2 = clf2.fit(mutTrainingData,resisTrainingData)

trainDatPredictions = list(clf2.predict(mutTrainingData))
testDatPredictions = list(clf2.predict(mutTestingData))

print(find_accuracy(trainDatPredictions, testDatPredictions))



def make_meshgrid(x, y, h=.02):
    """Create a mesh of points to plot in

    Parameters
    ----------
    x: data to base x-axis meshgrid on
    y: data to base y-axis meshgrid on
    h: stepsize for meshgrid, optional

    Returns
    -------
    xx, yy : ndarray
    """
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    return xx, yy


def plot_contours(ax, clf, xx, yy, **params):
    """Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)
    return out



X = mutationsDataFrame
Y = resistanceDataFrame
	
#Convert data frames to integer form
print("Converting data to integer form")
for i in range(len(X)):
     X[i] = list(map(float, X[i]))

for i in range(len(y)):
     y[i] = int(y[i][0])

X=np.array(X)
y = np.array(y)

# we create an instance of SVM and fit out data. We do not scale our
# data since we want to plot the support vectors

print("Building models, second round")
  # SVM regularization parameter
models = (svm.SVC(kernel='rbf', C=1),
          svm.SVC(kernel='poly', C=1, degree = 7))
models = (clf.fit(X, y) for clf in models)

print("Making Figure")
# title for the plots
titles = ('SVC with RBF kernel',
          'SVC with Polynomial kernel')

# Set-up 2x2 grid for plotting.
fig, sub = plt.subplots(2)
plt.subplots_adjust(hspace=0.4)

X0, X1 = X[:, 0], X[:, 1]
xx, yy = make_meshgrid(X0, X1)

for clf, title, ax in zip(models, titles, sub.flatten()):
    plot_contours(ax, clf, xx, yy,
                  cmap=plt.cm.coolwarm, alpha=0.8)
    ax.scatter(X0, X1, c=y, cmap=plt.cm.coolwarm, s=20, edgecolors='k')
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(title)
    ax.set_xticks(())
    ax.set_yticks(())

print("Showing plot")
plt.show()


