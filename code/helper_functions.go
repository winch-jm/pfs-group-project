package main

import (
	// "fmt"
)

// Test CSR Graphs
func DummyTriangle() *CSR {
    N := int32(3)

    // adjacency:
    // 0: 1,2
    // 1: 0,2
    // 2: 0,1

    indptr := []Idx{0, 2, 4, 6}          // len = N+1
    indices := []NodeID{1,2, 0,2, 0,1}   // neighbors concatenated
    data    := []Weight{1,1, 1,1, 1,1}

    degree := make([]float32, N)
    for i := int32(0); i < N; i++ {
        sum := float32(0)
        for k := indptr[i]; k < indptr[i+1]; k++ {
            sum += data[k]
        }
        degree[i] = sum
    }

    return &CSR{N: N, Indptr: indptr, Indices: indices, Data: data, Degree: degree}
}

func DummyTwoCliques() *CSR {
    N := int32(6)

    // adjacency: each clique is complete with weight=1
    // bridge edge 2—3 has weight=0.1
    
    indptr := []Idx{
        0,          // start row 0
        2,          // row 1
        4,          // row 2
        7,          // row 3
        10,          // row 4
        12,         // row 5
        14,         // end
    }

    indices := []NodeID{
        // row 0:
        1,2,
        // row 1:
        0,2,
        // row 2:
        0,1,3,
        // row 3:
        2,4,5,
        // row 4:
        3,5,
        // row 5:
        3,4,
    }

    data := []Weight{
        // clique A weights
        1,1, 
		1,1,  
		1,1,0.1,		// weak edge 2—3
        0.1,1,1,       // weak edge 3-2
        1,1, 
		1,1, 
    }

    degree := make([]float32, N)
    for i := int32(0); i < N; i++ {
        sum := float32(0)
		// fmt.Println(indptr[i],indptr[i+1])
        for k := indptr[i]; k < indptr[i+1]; k++ {
            sum += data[k]
        }
		// fmt.Println(sum,i)
        degree[i] = sum
    }

    return &CSR{N: N, Indptr: indptr, Indices: indices, Data: data, Degree: degree}
}