package main

import (
	"fmt"
  "math"
  "strconv"
  "io"
  "encoding/csv"
  "os"
  "log"
)


func CountCommunities(P Partition) int {
    seen := make(map[int32]struct{})
    for _, c := range P {
        seen[c] = struct{}{}
    }
    return len(seen)
}

func CSVToDenseRows(filename string) ([]string,DenseRows) {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Error opening file %q: %v", filename, err)
	}
	defer f.Close()

	reader := csv.NewReader(f)

	var dr DenseRows
	sigIds := make([]string,0)
	rowNum := 0
	for {
		// if rowNum > 500 {
		// 	break
		// }

		record, err := reader.Read()
		if rowNum > 0{
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("Error reading CSV: %v", err)
		}
		dr.D = len(record) - 2
		for ind , str := range record {
			if ind == 1 {
				sigIds = append(sigIds,str)

			} else if ind > 1{
			tempFloat, err := strconv.ParseFloat(str, 32)
			if err != nil {
				log.Fatalf("Error parsing %q as float32: %v", str, err)
			}
			dr.Data = append(dr.Data, float32(tempFloat))
		}
		}
	}
		rowNum ++
		dr.N++ // count rows
	}
	dr.N--
	fmt.Printf("Loaded matrix with %d rows and %d columns\n", dr.N, dr.D)
	fmt.Printf("Data slice length: %d (should be N*D = %d)\n",
		len(dr.Data), dr.N*dr.D)

	return sigIds, dr
}

func NormalizeDenseRows(dr DenseRows) DenseRows {
	fmt.Printf("Normalizing data (zero mean, dividing by L2 norm)...\n")
    // 1. Compute mean of each gene (column)
    means := make([]float32, dr.D)
    for gene := 0; gene < dr.D; gene++ {
        for idx := gene; idx < len(dr.Data); idx += dr.D {
            means[gene] += dr.Data[idx]
        }
        means[gene] /= float32(dr.N)
    }

    // 2. Subtract mean from each gene (column-center)
    for gene := 0; gene < dr.D; gene++ {
        for idx := gene; idx < len(dr.Data); idx += dr.D {
            dr.Data[idx] -= means[gene]
        }
    }

    // 3. L2-normalize each cell (row)
    for cell := 0; cell < dr.N; cell++ {
        start := cell * dr.D
        end := start + dr.D

        temp := dr.Data[start:end]
        diff := l2Norm(temp)
        if diff == 0 {
            continue // avoid division by zero
        }

        for i := start; i < end; i++ {
            dr.Data[i] /= diff
        }
    }

    return dr
}



func l2Norm(vector []float32) float32 {
	var sum float64
	for _, v := range vector {
		f := float64(v)
		sum += f * f
	}
	return float32(math.Sqrt(sum))
}
