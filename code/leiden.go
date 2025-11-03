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
	// return final Partition, stats
}

