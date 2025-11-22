package main

import (
	"container/heap"
	"crypto/x509"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
)

type DenseRows struct {
	// n rows (cells) × d columns (genes or features)
	N, D int
	Data []float32 // row-major: row i at Data[i*D:(i+1)*D]
}




// MaxHeap is a max-heap of ints.
type MaxHeap []int

func (h MaxHeap) Len() int           { return len(h) }
func (h MaxHeap) Less(i, j int) bool { return h[i] > h[j] } // NOTE: '>' for max-heap
func (h MaxHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }

// Push and Pop use pointer receiver because they modify the slice.
func (h *MaxHeap) Push(x any) {
	*h = append(*h, x.(int))
}

func (h *MaxHeap) Pop() any {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[:n-1]
	return x
}



func Cosine(a, b []float32) float32 {

	var s float32
	// dot product
	for i := range a {
		s += a[i]*b[i]
	}

	return s
}

func main() {
	file, err := os.Open("/Users/nethanramachandran/go/src/pfs-group-project/newCode/ctl_subset.csv")
	if err != nil {
		log.Fatalf("Error opening file: %v", err)
	}

	var file_dr DenseRows

	defer file.Close()
	reader := csv.NewReader(file)
	firstRow := true
	ctr := 0
	for {
		record, err := reader.Read()
		if ctr > 0 {

			if err == io.EOF {
				fmt.Println("There was an error reading the cell")
				break
			}

			if err != nil {
				log.Fatalf("Error reading CSV: %v", err)
			}

			if firstRow {
				file_dr.D = len(record)
				firstRow = false
			}

			for i := range record {
				str := record[i]
				tempFloat, err := strconv.ParseFloat(str, 32)
				fmt.Println(str)
				if err != nil {
					fmt.Println("there was an error parsing the string into a float32")
					break

				}
				file_dr.Data = append(file_dr.Data, float32(tempFloat))
			}
		}
		ctr++
		file_dr.N++
	}

	fmt.Printf("Loaded matrix with %d rows and %d columns\n", file_dr.N, file_dr.D)
	fmt.Printf("Data slice length: %d (should be N*D = %d)\n",
	len(file_dr.Data), file_dr.N*file_dr.D)





}

func Similarities(dr DenseRows) {
	


	for i:=0;i<dr.D;i++{
		
		for j:=0;j<dr.D;j++{
			tempCosineSimilarity:=make([]float32,0)
			make newMaxHeap(make)
			if i!=j{
				dataSet1:= dr.Data[i*dr.D:(i+1)*dr.D]
				dataSet2:= dr.Data[j*dr.D:(j+1)*dr.D]
				metric:= Cosine(dataSet1,dataSet2)
				tempCosineSimilarity = append(tempCosineSimilarity, metric)
			}
			// top3 he

			
		}
	}
}

func GraphCreation(){

}


// First thing: Turn CSV into denserows object
