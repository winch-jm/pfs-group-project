package main

import (
	// "fmt"
  "math"
)


func CountCommunities(P Partition) int {
    seen := make(map[int32]struct{})
    for _, c := range P {
        seen[c] = struct{}{}
    }
    return len(seen)
}

func CSVToDenseRows(file string){
  file, err := os.Open(file)
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
func normalizeDenseRows(dr DenseRows) DenseRows{
  means = make([]float32, dr.D)
  for gene := 0; gene < len(dr.D); gene++{
    for in := gene; in < len(dr.Data); in += dr.D{
      means[gene] = means[gene] + data[in]
    }
    means[gene] = means[gene] / dr.N
  }
  for gene := 0; gene < len(dr.D); gene++{
    for in := gene; in < len(dr.Data); in += dr.D{
      dr.Data[in] = dr.Data[in] - means[gene]
    }
  }
  for cell := 0; cell < len(dr.N); cell++{
    temp := dr.Data[cell * dr.D : (cell + 1) * dr.D]
    diff := l2Norm(temp)
    for ind := cell * dr.D; ind < (cell + 1) * dr.D; ind++{
      dr.Data[ind] = dr.Data[ind] / diff
    }
  }
return dr
}
func l2Norm(vector []float32) float32{
  sum := 0
  for x:= 0; x < len(vector); x++{
    sum = sum + (vector[x] * vector[x])
  }
  return math.Sqrt(sum)
}
