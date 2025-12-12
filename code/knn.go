// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

import (
	"container/heap"
   "fmt"
   "sync"
//    "runtime"
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
    fmt.Printf("Running Weight KNN:\n")
    fmt.Printf("Creating Adjacency Matrix..\n")
    // numWorkers := runtime.NumCPU()
	adjacencyMatrix := Similarities(dr, k)
    fmt.Printf("Constructing Graph...\n")
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

// RowResult holds the top-k neighbors for a given row i
type RowResult struct {
    i    int
    js   []int
    sims []float32
}

func SimilaritiesParallel(dr DenseRows, k int, numWorkers int) [][]float32 {
    n := dr.N

    // adjacency[i][j] = similarity between i and j (0 if not among top-k)
    adjacencyMatrix := make([][]float32, n)
    for i := range adjacencyMatrix {
        adjacencyMatrix[i] = make([]float32, n)
    }

    jobs := make(chan int, numWorkers)
    results := make(chan RowResult, numWorkers)

    var wg sync.WaitGroup

    // Worker goroutines
    worker := func() {
        defer wg.Done()

        for i := range jobs {
            rowI := dr.Data[i*dr.D : (i+1)*dr.D]

            scores := make([]*IndexScore, 0, n-1)
            for j := 0; j < n; j++ {
                if i == j {
                    continue
                }
                rowJ := dr.Data[j*dr.D : (j+1)*dr.D]
                sim := Cosine(rowI, rowJ)

                scores = append(scores, &IndexScore{
                    similarity: sim,
                    index:      int32(j),
                })
            }

            // Build max-heap
            h := MaxHeap(scores)
            heap.Init(&h)

            // Extract top-k
            js := make([]int, 0, k)
            sims := make([]float32, 0, k)
            for x := 0; x < k && h.Len() > 0; x++ {
                item := heap.Pop(&h).(*IndexScore)
                js = append(js, int(item.index))
                sims = append(sims, item.similarity)
            }

            // Send result for row i
            results <- RowResult{
                i:    i,
                js:   js,
                sims: sims,
            }
        }
    }

    // Start workers
    if numWorkers < 1 {
        numWorkers = 1
    }
    wg.Add(numWorkers)
    for w := 0; w < numWorkers; w++ {
        go worker()
    }

    // Feed jobs
    go func() {
        for i := 0; i < n; i++ {
            jobs <- i
        }
        close(jobs)
    }()

    // Close results once all workers are done
    go func() {
        wg.Wait()
        close(results)
    }()

    // Collect results and fill adjacency matrix (safely, in one goroutine)
    for res := range results {
        i := res.i
        for idx, j := range res.js {
            sim := res.sims[idx]
            adjacencyMatrix[i][j] = sim
            adjacencyMatrix[j][i] = sim // keep symmetric, same behavior as before
        }
    }

    return adjacencyMatrix
}
