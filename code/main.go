package main

import(
	"fmt"
	"log"
)

func main(){

	// Preprocess 
	// // - load dense rows
	file := "../data/MCF7_subset.csv"
	k := 25
	
	sigIds, dataset := CSVToDenseRows(file)

	// - L2 normalize
	dataset = NormalizeDenseRows(dataset)

	// build KNN (output is a CSR)
	g := WeightedKNN(dataset, k)

	// g := DummyTwoCliques()
	//free memory from denseRows (no longer needed)
	dataset = DenseRows{}
	//graph does not change after KNN!!

	// L
	P, levelStats := LeidenCommunityDetection(g, RBPM, 0.001, 1.0, 10, 10, 32)
	fmt.Println(len(P),P)
	fmt.Println(levelStats)

	g.WriteEdgeListCSV("../data/graph_edges_mcf7.csv",true)
	// ... run leiden, get origP ...
	if err := WriteFinalPartitionCSV(sigIds, P, "../data/partition_mcf7.csv"); err != nil {
		log.Fatal(err)
	}

}