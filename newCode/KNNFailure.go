
package main

import (
    "fmt"

    "github.com/sjwhitworth/golearn/base"
    "github.com/sjwhitworth/golearn/evaluation"
    "github.com/sjwhitworth/golearn/knn"
)

func main() {
    // Load CSV where the last column is the class label
    data, err := base.ParseCSVToInstances("data.csv", true)
    if err != nil {
        panic(err)
    }

    // Train/test split
    train, test := base.InstancesTrainTestSplit(data, 0.7)

    // New KNN classifier: k=5, use Euclidean distance
    cls := knn.NewKnnClassifier("euclidean", "linear", 5)
    // "linear" is the weighting scheme; can also use "inverse" or "gaussian"
    // depending on the GoLearn version you have

    cls.Fit(train)

    preds, err := cls.Predict(test)
    if err != nil {
        panic(err)
    }

    confMat, err := evaluation.GetConfusionMatrix(test, preds)
    if err != nil {
        panic(err)
    }

    fmt.Println(evaluation.GetSummary(confMat))
}