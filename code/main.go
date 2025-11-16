package main

import(
	"fmt"
)

func main(){
	// Test Local Moving and Refinement

	g := DummyTwoCliques()
	cs := InitializeCommStats(g)
	mb := InitializeMoveBuffers(32)
	P := InitializePartition(g)

	P_initial := make(Partition,int(g.N))
	copy(P_initial,P)


	cfg := QualityCfg{Fn:RBPM,SelfLoop:0.0,Gamma:0.5}

	// Q0 := ComputeModularity(g, P)
	_,_ = LocalMove(g,cfg,P,cs,mb,10)
	P_unrefined := P
	rb := InitializeRefineBuffers(g,32)

	P_refined := RefinePartition(g, P, rb)

	fmt.Println("Initial Partition: ", P_initial)
	fmt.Println("After Local Moving: ", P_unrefined)
	fmt.Println("After Refinement: ", P_refined)


}