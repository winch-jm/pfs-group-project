// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

import(
	"fmt"
	"log"
	"strconv"
	"os"
)

func RunGridSearch(dataset DenseRows, saveDir, prefix string, kList []int, resList []float32, sigIds []string) {
	
	for _,k := range kList {

		// build KNN (output is a CSR)
		g := WeightedKNN(dataset, k)

		for _, gamma := range resList {

			// Leiden algorithm
			P, levelStats := LeidenCommunityDetection(g, RBPM, 0.001, gamma, 10, 10, 32)
			
			fmt.Println(levelStats)

			// Write partition to file
			if err := WriteFinalPartitionCSV(sigIds, P, saveDir+"/gridsearch/"+prefix+"_partition_k" + strconv.Itoa(k)+"_g"+strconv.Itoa(int(gamma*100))+".csv"); err != nil {
				log.Fatal(err)
			}
		}
		// Write graph edges to file
		g.WriteEdgeListCSV(saveDir +"/gridsearch/"+prefix+"_graph_edges_k"+ strconv.Itoa(k)+".csv",true)
	}
}

func main() {
	//file := "../data/mcf7_subset_pca.csv"
	file := os.Args[1] //First argument is an input of CSV file
	projectName := os.Args[2]//Second argument is the prefix tagged onto images
	outputDir := os.Args[3]//Third argument is the output directory for graphs/partition labels
	gridSearch := os.Args[4]//Fourth argument is whether user wants to run GridSearch to find optimal k and gamma

	// read in dataset
	sigIds, dataset := CSVToDenseRows(file)

	// normalize signatures/pca vectors
	dataset = NormalizeDenseRows(dataset)

	// run grid search optimization
	if gridSearch == "true"{ //Runs grid search and then runs the algorithm
		err := os.MkdirAll(outputDir+"/gridsearch/", 0755) 
		if err != nil {
			fmt.Printf("Error creating directory: %s\n", err)
			return
		}
		kList := []int{10,15,20,25} 
		resList := []float32{0.25,0.5,0.75,1.0,1.25,1.5}
		RunGridSearch(dataset, outputDir, projectName, kList, resList, sigIds)

	} else { //If they don't specify gridSearch, user has to define k and gamma and run algorithm

		err := os.MkdirAll(outputDir, 0755) 
		if err != nil {
			fmt.Printf("Error creating directory: %s\n", err)
			return
		}

		k, _ := strconv.Atoi(os.Args[5]) //
		gamma, _ := strconv.ParseFloat(os.Args[6], 64)
		g := WeightedKNN(dataset, k)
		P, levelStats := LeidenCommunityDetection(g, RBPM, 0.001, float32(gamma), 10, 10, 32)
		fmt.Println(levelStats)

		// Write partition to file
		if err := WriteFinalPartitionCSV(sigIds, P, outputDir+"/gridsearch/"+projectName+"_partition_k" + strconv.Itoa(k)+"_g"+strconv.Itoa(int(gamma*100))+".csv"); err != nil {
			log.Fatal(err)
		}
		
		// Write graph edges to file
		g.WriteEdgeListCSV(outputDir +"/gridsearch/"+projectName+"_graph_edges_k"+ strconv.Itoa(k)+".csv",true)
	}
}

//command line parsing for input dataset,output directory if you are doing grid_search, if you want to run algorithm once specify k and gamma
//