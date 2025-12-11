// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

import(
	// "fmt"
	// "log"
	// "strconv"
)

func main(){

	// Preprocess 
	// // - load dense rows
	file := "../data/mcf7_subset_pca.csv"

	kList := []int{10,15,20,25} 
	// resList := []float32{0.25,0.5,0.75,1.0,1.25,1.5}
	for _,k := range kList {
		_, dataset := CSVToDenseRows(file)
			
		// - L2 normalize
		dataset = NormalizeDenseRows(dataset)

		// build KNN (output is a CSR)

		WeightedKNN(dataset, k)
		// fmt.Println(g)
		
	// 	//free memory from denseRows (no longer needed)
	// 	dataset = DenseRows{}
	// 	//graph does not change after KNN!!

	// 	for _, gamma := range resList {

	// 		// Leiden algorithm
	// 		P, levelStats := LeidenCommunityDetection(g, RBPM, 0.001, gamma, 10, 10, 32)
	// 		fmt.Println(len(P),P)
	// 		fmt.Println(levelStats)

			
	// 		// ... run leiden, get origP ...
	// 		if err := WriteFinalPartitionCSV(sigIds, P, "../data/gridsearch/partition_mcf7_k" + strconv.Itoa(k)+"_g"+strconv.Itoa(int(gamma*100))+".csv"); err != nil {
	// 			log.Fatal(err)
	// 		}
			
	// 	}
	// 	g.WriteEdgeListCSV("../data/gridsearch/graph_edges_mcf7_k"+ strconv.Itoa(k)+".csv",true)
	}
	

}