package main

import (
	"math"
)

// cosine simlarity between two normalized vectors
func Cosine(a, b []float32) float32 {

	var s float32
	// dot product
	for i := range a {
		s += a[i]*b[i]
	}

	return s
}



// Leiden Metrics

func Modularity(A Matrix) float64 {
	m := SumEdgeWeights()
	s := 0
	for i := range A {
		for j := range A {
			s += (A[i][j] - (sum_i*sum_j)/(2*m))*KroneckerDelta(ci,cj)
		}
	}

	return 1.0/(2.0*m) * s 

}

// Kronecker Delta function
func KroneckerDelta(c1, c1 *Node) {
	if c1.community == c2.community {
		return 1.0
	}
	return 0.0
}

// Reichardt Bornholdt Potts Model
func RBPM(A Matrix, gamma float64) float64 {
	m := SumEdgeWeights()
	s := 0
	for i := range A {
		for j := range A {
			s += (A[i][j] - gamma*(sum_i*sum_j)/(2*m))*KroneckerDelta(ci,cj)
		}
	}

	return s 

}

// Constant Potts Model
func CPM(A Matrix, gamma float64) float64 {
	m := SumEdgeWeights()
	s := 0
	for i := range A {
		for j := range A {
			s += (A[i][j] - gamma)*KroneckerDelta(ci,cj)
		}
	}

	return s 
}
