package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
)

// Example represents one row of the dataset.
type Example struct {
	Features []float64
	Label    string
}

// ----------------- Distance & weights -----------------

// Euclidean distance between two feature vectors.
func euclidean(a, b []float64) float64 {
	if len(a) != len(b) {
		log.Fatalf("dimension mismatch: %d vs %d", len(a), len(b))
	}
	var sum float64
	for i := range a {
		d := a[i] - b[i]
		sum += d * d
	}
	return math.Sqrt(sum)
}

// Inverse-distance weight (with epsilon to avoid division by zero).
func inverseDistance(dist float64) float64 {
	const eps = 1e-9
	if dist < eps {
		dist = eps
	}
	return 1.0 / dist
}

// ----------------- kNN prediction -----------------

// predictOne predicts the label for example at index idx, using all other
// examples in train as neighbors, with weighted kNN.
func predictOne(train []Example, idx int, k int) string {
	query := train[idx].Features

	// Collect distances to all other points.
	type neighbor struct {
		dist  float64
		label string
	}
	nbs := make([]neighbor, 0, len(train)-1)
	for j, ex := range train {
		if j == idx {
			continue // don't use the point itself as a neighbor
		}
		d := euclidean(query, ex.Features)
		nbs = append(nbs, neighbor{dist: d, label: ex.Label})
	}

	// Select k nearest neighbors (simple O(n*k) partial selection).
	if k > len(nbs) {
		k = len(nbs)
	}
	for i := 0; i < k; i++ {
		best := i
		for j := i + 1; j < len(nbs); j++ {
			if nbs[j].dist < nbs[best].dist {
				best = j
			}
		}
		nbs[i], nbs[best] = nbs[best], nbs[i]
	}

	// Weighted vote by inverse distance.
	score := make(map[string]float64)
	for i := 0; i < k; i++ {
		w := inverseDistance(nbs[i].dist)
		score[nbs[i].label] += w
	}

	// Pick label with highest weight.
	var bestLabel string
	bestScore := math.Inf(-1)
	for lab, s := range score {
		if s > bestScore {
			bestScore = s
			bestLabel = lab
		}
	}
	return bestLabel
}

// ----------------- CSV loading -----------------

// loadCSV loads a CSV where the last column is the label, the others are numeric features.
func loadCSV(path string) ([]Example, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := csv.NewReader(f)
	records, err := r.ReadAll()
	if err != nil {
		return nil, err
	}
	if len(records) < 2 {
		return nil, fmt.Errorf("need at least header + 1 row")
	}

	// Skip header, assume last column is label.
	var data []Example
	for i, row := range records[1:] {
		if len(row) < 2 {
			return nil, fmt.Errorf("row %d has fewer than 2 columns", i+2)
		}
		features := make([]float64, len(row)-1)
		for j := 0; j < len(row)-1; j++ {
			v, err := strconv.ParseFloat(row[j], 64)
			if err != nil {
				return nil, fmt.Errorf("row %d col %d not numeric: %v", i+2, j+1, err)
			}
			features[j] = v
		}
		label := row[len(row)-1]
		data = append(data, Example{Features: features, Label: label})
	}
	return data, nil
}
