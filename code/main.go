package main

import(
	"fmt"
)

func main(){

	// // Preprocess 
	// // - load dense rows
	// dataset := InitializeDenseRows()
	// // - L2 normalize
	// dataset := L2Normalize(dataset)
	// // build KNN (output is a CSR)
	// g := WeightedKNN(dataset)
	// free memory from denseRows (no longer needed)
	// graph does not change after KNN!!

	// Test Local Moving and Refinement

	g := DummyTwoCliques()
	modTypes := []QualityFn{Modularity,RBPM,CPM}

	for _,modFn := range modTypes {
		P, levelStats := LeidenCommunityDetection(g, modFn, 0.0, 0.5, 10, 10, 32)
		fmt.Println(P)
		fmt.Println(levelStats)
	}
	
	
	// P_initial := make(Partition,int(g.N))
	// copy(P_initial,P)


	

	// // Q0 := ComputeModularity(g, P)
	// _,_ = LocalMove(g,cfg,P,cs,mb,10)
	// P_unrefined := P
	

	// P_refined := RefinePartition(g, P, rb)

	// fmt.Println("Initial Partition: ", P_initial)
	// fmt.Println("After Local Moving: ", P_unrefined)
	// fmt.Println("After Refinement: ", P_refined)


}