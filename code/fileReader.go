package main

import (
	"encoding/csv"
	"log"
	"os"
	"strconv"
)

func ReadCSV(file string) []*Cell {
	// creates a list of pointers to cells
	cells := make([]*Cell, 0)
	// opens the file
	f, err := os.Open("data/ctl_subset.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	// creates a reader
	reader := csv.NewReader(f)
	header, err := reader.Read()
	records, err := reader.ReadAll()
	// creates a map with a key as a string and a value as a float64
	// is it neccesary to create it twice
	temp_dict := make(map[string]float64)
	// ranges through each row
	for i, row := range records {
		temp_dict = make(map[string]float64)
		// never takes in the first row or column because those are headers
		if i > 0 {
			for ind := range row {
				if ind != 0 {
					// if the row and column > 0, then make the key the header index,
					// and the value the parsed row index.
					temp_dict[header[ind]], _ = strconv.ParseFloat(row[ind], 64)
				}
			}
			// once the map has been c
			cells = append(cells, &Cell{row[0], temp_dict})
		}
	}
	return cells
}
