package main

// https://en.wikipedia.org/wiki/Leiden_algorithm
import (
	"fmt"
)


// Leiden skeleton
// func LeidenCommunityDetection(dataset string) (Partition, LevelStats) {
// 	// Preprocess 
// 	// - load dense rows
// 	dataset := InitializeDenseRows()
// 	// - L2 normalize
// 	dataset := L2Normalize(dataset)
// 	// build KNN (output is a CSR)
	
// 	G := WeightedKNN(dataset)
// 	// free memoryd from denseRows (no longer needed)
// 	// graph does not change after KNN!!


// 	// Initialize
// 	// - initialize partition
// 	P := InitializePartition(dataset)
// 	// - initialize community statistics
// 	cs := InitializeCommStats()
// 	// - initiialize movebuffers
// 	mb := InitializeMoveBuf()
// 	// - initiialize refinebuffers
// 	rb := InitializeRefineBuf()



// 	// Iterate for maxLevels
// 	for level := 0; level < maxLevels; level++ {
// 		// 1) Local moving (optimize quality with node moves)
// 		moves, deltaQ = LocalMove()
//		// termination condition
// 		if moves == 0 && math.Abs(deltaQ) < eps {
// 			break
// 		}
// 		// 2) Refinement (split poorly connected communities)
// 		Refine()

// 		// 3) Aggregation (contract communities --> supernodes)
// 		Aggregate()

// 		// Bookkeeping

// 		
// 	}
// 	// return final Partition, stats
// }

// LocalMove runs the Leiden/Louvain local moving phase.
// Returns total accepted moves and total DeltaQ gain across all sweeps.
func LocalMove(g *CSR, cfg QualityCfg, P Partition, cs *CommStats, mb *MoveBuffers, maxSweeps int) (int, float32) {
	if maxSweeps <= 0 { maxSweeps = 1<<30 } // effectively "until convergence"

	totalMoves := 0
	totalGain  := float32(0.0)

	// Node order: deterministic (0..N-1). For randomness, shuffle an index slice each sweep.
	order := make([]int32, int(g.N))
	for i := range order { order[i] = int32(i) }

	for sweep := 0; sweep < maxSweeps; sweep++ {
		// get # of moves and deltaQ gains
		movesThis, gainsThis := LocalMoveSweep(g, cfg, P, cs, mb, order)

		totalMoves += movesThis
		totalGain  += gainsThis

		if movesThis == 0 { break } // converged
	}
	
	return totalMoves, totalGain
}

// Run local move for one pass (sweep) over all nodes
func LocalMoveSweep(
	g *CSR,
	cfg QualityCfg,
	P Partition,
	cs *CommStats,
	mb *MoveBuffers,
	order []int32) (int, float32) {

	movesThis, gainsThis := 0, float32(0.0)
	
	for _, i := range order {
		currComm := P[i]
		ki := g.Degree[i]

		// remove node i from current community
		cs.EnsureCapacity(currComm)

		cs.Tot[currComm] -= ki
		cs.Size[currComm]--

		
		// accumulate k_i,in(c) over neighbor communities
		mb.Reset()
		for idx := g.Indptr[i]; idx < g.Indptr[i+1]; idx++ {
			j := g.Indices[idx]
			w := g.Data[idx]
			mb.Add(P[j], w)
		}

		// self-loop into "stay" community comm
		if cfg.SelfLoop > 0 {
			if idx, ok := mb.Seen[currComm]; ok {
				mb.CommWts[idx] += float32(cfg.SelfLoop)
			} else {
				mb.Add(currComm, float32(cfg.SelfLoop))
			}	
		}

		// Evaluate ΔQ for each candidate community c
		bestC := currComm
		bestDelta := float32(0.0)

		for t:= 0; t < len(mb.CommIDs); t++ {
			c := mb.CommIDs[t]
			kin := float32(mb.CommWts[t])

			cs.EnsureCapacity(c)

			var deltaQ float32
			switch cfg.Fn {
			case Modularity:
				deltaQ = deltaQModularity(kin, ki, cs.Tot[c], g.TwoM)
			case RBPM:
				deltaQ = deltaRBPM(kin, ki, cs.Tot[c], cfg.Gamma, g.TwoM)
			case CPM:
				deltaQ = deltaCPM(kin, cfg.Gamma, cs.Size[c])
			default:
				deltaQ = deltaQModularity(kin, ki, cs.Tot[c], g.TwoM)
			}

			if deltaQ > bestDelta || (deltaQ == bestDelta && c < bestC) {
				bestDelta = deltaQ
				bestC = c
			}

		}

		Qbefore := ComputeModularity(g, P)  // uses float64, full recompute

		// Commit: move if positive gain; else reinsert into comm
		if bestDelta > 0 && bestC != currComm {
			P[i] = bestC
			cs.Tot[bestC] += float32(ki)
			cs.Size[bestC]++

			Qafter := ComputeModularity(g, P)

			fmt.Printf("Move node %d: %d -> %d, bestDelta=%.6f, Qdiff=%.6f\n",
				i, currComm, bestC, bestDelta, Qafter-Qbefore)
			movesThis++
			gainsThis += bestDelta
		} else {
			// revert to currComm
			cs.Tot[currComm] += float32(ki)
			cs.Size[currComm]++ 
		}

	}
	return movesThis, gainsThis
}

// Split communities if they contain multiple connected components, relabel Partition with consecutive labels starting at 0
func RefinePartition(g *CSR, P Partition, rb *RefineBuffers) Partition {
	n := int(g.N)

	refined := make(Partition, n)

	for i := range refined {
		refined[i] = -1 // marks as "unassigned"
	}

	var nextComm int32

	for start := int32(0); start < g.N; start++ {
		if refined[start] != -1 {
			continue // already assigned to some refined community
		}


		origComm := P[start]

		// start new BFS for connected component
		rb.Reset()	// clears queue, increments epoch counter
		rb.Queue = append(rb.Queue,start)
		rb.Visited[start] = rb.Stamp
		refined[start] = nextComm

		for len(rb.Queue) > 0 {
			// pop first element
			v := rb.Queue[0]
			rb.Queue = rb.Queue[1:]

			rowStart := g.Indptr[v]
			rowEnd := g.Indptr[v+1]

			for idx := rowStart; idx < rowEnd; idx++ {
				u := g.Indices[idx]

				// only traverse edges within same original community
				if P[u] != origComm {
					continue
				}

				// skip if already assigned to this or another refined community
				if refined[u] != -1 {
					continue
				}

				// already visited within BFS?
				if rb.Visited[u] == rb.Stamp {
					continue
				}

				rb.Visited[u] = rb.Stamp
				refined[u] = nextComm
				rb.Queue = append(rb.Queue,u)
			}

		}
		nextComm++
	}

	return refined
}

// func Aggegrate() {
	
// }