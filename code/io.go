package main

import (
	"encoding/csv"
	"log"
	"os"
	"strconv"
)

func ReadCSV(file string) []*Cell{
	cells := make([]*Cell, 0)
	f, err := os.Open("data/ctl_subset.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	reader := csv.NewReader(f)
	header, err := reader.Read()
	records, err := reader.ReadAll()
	temp_dict := make(map[string]float64)
	for i, row := range records {
		temp_dict = make(map[string]float64)
		if i > 0{
		for ind, _ := range row{
			if ind != 0{
			temp_dict[header[ind]], _ = strconv.ParseFloat(row[ind], 32)
			}
		}
		cells = append(cells,&Cell{row[0], temp_dict})
		}
	}
	return cells
}