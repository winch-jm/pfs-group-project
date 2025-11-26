package main

import(
	"fmt"
	//"log"
)

func main(){

	// Preprocess 
	// - load dense rows
	file := "ctl_subset.csv"
	k := 30
	// dataset := CSVToDenseRows(file)
	dataset := CSVToDenseRows(file)
	fmt.Println("hello")
	// - L2 normalize
	dataset = NormalizeDenseRows(dataset)
	// // build KNN (output is a CSR)
	g := WeightedKNN(dataset, k)

	//free memory from denseRows (no longer needed)
	//graph does not change after KNN!!

	//Test Local Moving and Refinement
	// g := DummyTwoCliques()
	// pG := Partition{0,0,0,1,1,1}
	// Aggregation(g, pG)

	modTypes := []QualityFn{Modularity,RBPM,CPM}
	var P Partition
	for _,modFn := range modTypes {
		P, levelStats := LeidenCommunityDetection(g, modFn, 0.0, 0.5, 10, 10, 32)
		fmt.Println(P)
		fmt.Println(levelStats)
	}
	
	
	// P_initial := make(Partition,int(g.N))
	// copy(P_initial, P)


	

	// // Q0 := ComputeModularity(g, P)
	// _,_ = LocalMove(g,cfg,P,cs,mb,10)
	// P_unrefined := P
	

	// P_refined := RefinePartition(g, P, rb)

	g_new := Aggregation(g, P)
	fmt.Println(g_new)


	// fmt.Println("Initial Partition: ", P_initial)
	// fmt.Println("After Local Moving: ", P_unrefined)
	// fmt.Println("After Refinement: ", P_refined)
}