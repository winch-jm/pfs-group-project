package main

import (
	"container/heap"
   //"fmt"
	//"strconv"
)





// MaxHeap is a max-heap of ints.

type MaxHeap []*IndexScore

func (h MaxHeap) Len() int           { return len(h) }
func (h MaxHeap) Less(i, j int) bool { return h[i].similarity > h[j].similarity }
func (h MaxHeap) Swap(i, j int) {
    h[i], h[j] = h[j], h[i]
    h[i].heapIdx = i
    h[j].heapIdx = j
}

func (h *MaxHeap) Push(x any) {
    item := x.(*IndexScore)
    item.heapIdx = len(*h)
    *h = append(*h, item)
}

func (h *MaxHeap) Pop() any {
    old := *h
    n := len(old)
    item := old[n-1]
    *h = old[:n-1]
    item.heapIdx = -1 // optional: mark as "not in heap" 
    return item
}
func WeightedKNN(dr DenseRows, k int) *CSR{
	adjacencyMatrix := Similarities(dr, k)
	return GraphCreation(adjacencyMatrix, dr)
}

func Similarities(dr DenseRows, k int) [][]float32 {
    n := dr.N

    // adjacency[i][j] = similarity between i and j (0 if not among top-k)
    adjacencyMatrix := make([][]float32, n)
    for i := range adjacencyMatrix {
        adjacencyMatrix[i] = make([]float32, n)
    }

    for i := 0; i < n; i++ {
        scores := make([]*IndexScore, 0, n-1)

        rowI := dr.Data[i*dr.D : (i+1)*dr.D]
        for j := 0; j < n; j++ {
            if i == j {
                continue
            }
            rowJ := dr.Data[j*dr.D : (j+1)*dr.D]
            sim := Cosine(rowI, rowJ) // you should already have this function
            scores = append(scores, &IndexScore{
                similarity: sim,
                index:  int32(j),
            })
        }

        // Build max-heap
        h := MaxHeap(scores)
        heap.Init(&h)

        // Take top-k neighbors
        for x := 0; x < k && h.Len() > 0; x++ {
            item := heap.Pop(&h).(*IndexScore)
            j := item.index
            adjacencyMatrix[i][j] = item.similarity
            adjacencyMatrix[j][i] = item.similarity // make it symmetric
        }
    }

    return adjacencyMatrix
}

func GraphCreation(ad [][]float32, dr DenseRows) *CSR {
    n := dr.N

    newIndptr := make([]Idx, n+1)
    newIndices := make([]NodeID, 0)
    newData := make([]Weight, 0)
    newDegree := make([]float32, n)
    totalEdges := 0
    var twoM float32 // keep this as float32 to match CSR.TwoM

    newIndptr[0] = 0

    for cell := 0; cell < n; cell++ {
        row := ad[cell]
        var deg float32
        for j, w := range row {
            if w != 0 {
                newIndices = append(newIndices, NodeID(j))
                newData = append(newData, Weight(w))
                deg += w
                totalEdges++
            }
        }

        newIndptr[cell+1] = Idx(totalEdges)
        newDegree[cell] = deg
        twoM += deg
    }

    return &CSR{
        N:       int32(n),      // number of nodes = number of rows (cells)
        Indptr:  newIndptr,
        Indices: newIndices,
        Data:    newData,
        Degree:  newDegree,
        TwoM:    twoM,
    }
}
// func AddandCountNonZero(vals []float32, newIndices []NodeID, newData []Weight)(int, float32){
// 	sm := 0
// 	degree := 0.0
// 	for x := range vals{
// 		if vals[x] != 0{
// 			sm ++
// 			newIndices = append(newIndices, NodeID(x))
// 			newData = append(newData, Weight(vals[x]))
// 			degree = float32(degree) + vals[x]

// 		}
// 	}
// 	return sm, degree
// }