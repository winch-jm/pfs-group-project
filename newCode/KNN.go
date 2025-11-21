package main

import (
	// "container/heap"
	// "crypto/x509"
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

// First thing: Turn CSV into denserows object
