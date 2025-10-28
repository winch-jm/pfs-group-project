package main

// https://en.wikipedia.org/wiki/Leiden_algorithm

func LeidenCommunityDetection(A Matrix) {
	for len(P) != len(A) {
		P := FastLouvainMoveNodes(G, P)

		P_refined := GetPRefined(G,P)

		G := AggregateGraph(G, P_refined)
	}
	

}

func FastLouvainMoveNodes() {

}

func GetPRefined() {

}

func RefinePartitionSubset() {

}

func AggregateGraph() {

}