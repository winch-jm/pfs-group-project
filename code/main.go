// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

import(
	// "fmt"
	"log"
	"strconv"
	// "os"
)

func RunGridSearch(dataset DenseRows, saveDir, prefix string, kList []int, resList []float32, sigIds []string) {
	
	for _,k := range kList {

		// build KNN (output is a CSR)
		g := WeightedKNN(dataset, k)

		for _, gamma := range resList {

			// Leiden algorithm
			P, levelStats := LeidenCommunityDetection(g, RBPM, 0.001, gamma, 10, 10, 32)
			print(levelStats)
			if err := WriteFinalPartitionCSV(sigIds, P, saveDir+"/gridsearch/"+prefix+"_partition_k" + strconv.Itoa(k)+"_g"+strconv.Itoa(int(gamma*100))+".csv"); err != nil {
				log.Fatal(err)
			}
		}
		g.WriteEdgeListCSV(saveDir +"/gridsearch/"+prefix+"_graph_edges_k"+ strconv.Itoa(k)+".csv",true)
	}
}

func main() {
	// params := os.Args
	file := "../data/mcf7_subset_pca.csv"
	prefix := "mcf7"
	// ./code filename 
	// fmt.Println(params)
	
	sigIds, dataset := CSVToDenseRows(file)
			
		// - L2 normalize
	dataset = NormalizeDenseRows(dataset)
	// Preprocess 
	// // - load dense rows
	
	kList := []int{10,15,20,25} 
	resList := []float32{0.25,0.5,0.75,1.0,1.25,1.5}
	saveDir := ""
	RunGridSearch(dataset, saveDir, prefix, kList, resList, sigIds)

}