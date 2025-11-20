package main

// https://en.wikipedia.org/wiki/Leiden_algorithm
import (
	"math"
)

func Leiden() {
	// Preprocess 
	// - load dense rows
	// - L2 normalize
	// build KNN (output is a CSR)
	// free memory from denseRows (no longer needed)
	// graph does not change after KNN!!


	// Initialize
	// - initialize partition
	// - initialize community statistics
	// - initiialize movebuffers
	// - initiialize refinebuffers



	// Iterate for maxLevels
	for level := 0; level < maxLevels; level++ {
		// 1) Local moving (optimize quality with node moves)
		moves, deltaQ = LocalMove()

		// 2) Refinement (split poorly connected communities)

		// 3) Aggregation (contract communities --> supernodes)

		// Bookkeeping

		// termination condition
		if moves == 0 && math.Abs(deltaQ) < eps {
			break
		}
	}
	// return final Partition, statsi
}

func Aggregation(graph *CSR, partitionGraph *Partition){
	N := int(graph.N)
    if len(partitionGraph) != N {
        panic("Aggregation: partition length != number of nodes")
    }

    // --- 1. Get number of communities C from partitionGraph ---
    maxComm := int32(0)
    for _, cid := range partitionGraph {
        if cid > maxComm {
            maxComm = cid
        }
    }
    C := int(maxComm) + 1 // communities are 0..C-1
    newN := C

    // nodeToComm[i] is just partitionGraph[i]
    nodeToComm := partitionGraph

    // --- 2. Aggregate edge weights between communities ---
    type edgeKey struct {
        u, v int32
    }
    edgeWeights := make(map[edgeKey]Weight)

    for u := 0; u < N; u++ {
        cu := nodeToComm[u]
        rowStart := int(graph.Indptr[u])
        rowEnd := int(graph.Indptr[u+1])

        for e := rowStart; e < rowEnd; e++ {
            v := int(graph.Indices[e])
            if u > v { // skip symmetric copy (undirected CSR)
                continue
            }
            cv := nodeToComm[v]

            a, b := cu, cv
            if a > b {
                a, b = b, a
            }
            k := edgeKey{u: a, v: b}
            edgeWeights[k] += graph.Data[e]
        }
    }

    // --- 3. Convert aggregated map into adjacency lists ---
    type nbr struct {
        to int32
        w  Weight
    }
    adj := make([][]nbr, newN)

    for k, w := range edgeWeights {
        cu, cv := k.u, k.v
        if cu == cv {
            // self-loop
            adj[cu] = append(adj[cu], nbr{to: cv, w: w})
        } else {
            // undirected: store both directions
            adj[cu] = append(adj[cu], nbr{to: cv, w: w})
            adj[cv] = append(adj[cv], nbr{to: cu, w: w})
        }
    }

    // --- 4. Build CSR arrays for new graph ---
    newIndptr := make([]Idx, newN+1)
    totalEdges := 0
    for i := 0; i < newN; i++ {
        newIndptr[i] = Idx(totalEdges)
        totalEdges += len(adj[i])
    }
    newIndptr[newN] = Idx(totalEdges)

    newIndices := make([]NodeID, totalEdges)
    newData := make([]Weight, totalEdges)

    pos := 0
    for i := 0; i < newN; i++ {
        for _, e := range adj[i] {
            newIndices[pos] = NodeID(e.to)
            newData[pos] = e.w
            pos++
        }
    }

    // --- 5. Degrees for each super-node ---
    newDegree := make([]float32, newN)
    for i := 0; i < newN; i++ {
        sum := float32(0)
        start := int(newIndptr[i])
        end := int(newIndptr[i+1])
        for k := start; k < end; k++ {
            sum += float32(newData[k])
        }
        newDegree[i] = sum
    }

    return &CSR{
        N:      int32(newN),
        Indptr: newIndptr,
        Indices: newIndices,
        Data:   newData,
        Degree: newDegree,
    }
	//Use partition list to compress CSR
	//AggMap is essentially a snapshot of each level's paritionGraph, which will end up getting smaller

}
