package main

// https://en.wikipedia.org/wiki/Leiden_algorithm
import (
	"fmt"
	"math"
)


// Leiden Community Detection: https://www.nature.com/articles/s41598-019-41695-z
func LeidenCommunityDetection(g *CSR, qualityFn QualityFn, selfLoop, gamma float32, maxSweeps, maxLevels, capGuess int) (Partition, []LevelStats) {

	// Initialize
	// - initialize partition
	P := InitializePartition(g)
	// - initiialize movebuffers
	mb := InitializeMoveBuffers(capGuess)
	// - initiialize refinebuffers
	rb := InitializeRefineBuffers(g, capGuess)

	// modularity/resolution params
	cfg := QualityCfg{Fn:qualityFn,SelfLoop:selfLoop,Gamma:gamma}

	statsList := make([]LevelStats,0)

	
	// Iterate for maxLevels
	for level := 0; level < maxLevels; level++ {
		// reinitialize CommStats info every level
		cs := InitializeCommStatsFromPartition(g,P)

		// 1) Local moving (optimize quality with node moves)
		moves, deltaQ := LocalMove(g, cfg, P, cs, mb, maxSweeps)

		// termination condition
		if moves == 0 && math.Abs(float64(deltaQ)) < 1e-4 {
			break
		}

		// 2) Refinement (split poorly connected communities)
		P = RefinePartition(g, P, rb)

		// 3) Aggregation (contract communities --> supernodes)
		// Aggregate()

		// Bookkeeping
		levelStats := LevelStats{
			Level:level, 
			Quality:ComputeModularity(g,P), // need to implement other qualityFn for this
			NumCommunities:CountCommunities(P), 
			Moves:moves,
		}
		statsList = append(statsList,levelStats)

	}
	return P, statsList
}

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

func Aggregation(graph *CSR, partitionGraph Partition) *CSR {
    N := int(graph.N)
    if len(partitionGraph) != N {
        panic("Aggregation: partition length != number of nodes")
    }

    // --- 1. Get number of communities C from partitionGraph ---
    maxComm := int32(0)
    for _, cid := range partitionGraph {
        if cid > maxComm {
            maxComm = cid
        }
    }

    newN := int(maxComm) + 1 // communities are 0..C-1
    // nodeToComm[i] is just partitionGraph[i]
    nodeToComm := partitionGraph

    // --- 2. Aggregate edge weights between communities ---
    type edgeKey struct {
        u, v int32
		w float32
    }

	//defines edges between two different communities
    edgeWeights := make(map[edgeKey]Weight)

    for u := 0; u < N; u++ {
        cu := nodeToComm[u]
        rowStart := int(graph.Indptr[u])
        rowEnd := int(graph.Indptr[u+1])

        for e := rowStart; e < rowEnd; e++ {
            v := int(graph.Indices[e])
            if u > v { 
                continue
            }
            cv := nodeToComm[v]

            a, b := cu, cv
            if a > b {
                a, b = b, a
            }
            k := edgeKey{u: a, v: b}
            edgeWeights[k] += graph.Data[e]
        }
    }

    // --- 3. Convert aggregated map into adjacency lists ---
	//Fix this
	fmt.Println(edgeWeights)
    adj := make([][]edgeKey, newN)

    for k, w := range edgeWeights {
        cu, cv := k.u, k.v
		if cu != cv{
            adj[cu] = append(adj[cu], edgeKey{u: cu, v: cv, w: w})
            adj[cv] = append(adj[cv], edgeKey{u: cv, v: cu, w: w})
        }
    }
	fmt.Println(adj)
	fmt.Println(newN)
    // --- 4. Build CSR arrays for new graph ---
    newIndptr := make([]Idx, newN+1)
    totalEdges := 0
    for i := 0; i < newN; i++ {
        newIndptr[i] = Idx(totalEdges)
        totalEdges += len(adj[i])
    }
    newIndptr[newN] = Idx(totalEdges)

    newIndices := make([]NodeID, totalEdges)
    newData := make([]Weight, totalEdges)

    pos := 0
    for i := 0; i < newN; i++ {
        for _, edge := range adj[i] {
            newIndices[pos] = NodeID(edge.v)
            newData[pos] = edge.w
            pos++
        }
    }

    // --- 5. Degrees for each super-node ---
    newDegree := make([]float32, newN)
	var twoM float32 
    for i := 0; i < newN; i++ {
        sum := float32(0)
        start := int(newIndptr[i])
        end := int(newIndptr[i+1])
        for k := start; k < end; k++ {
            sum += float32(newData[k])
        }
        newDegree[i] = sum
		twoM += sum
    }

    return &CSR{
        N:      int32(newN),
        Indptr: newIndptr,
        Indices: newIndices,
        Data:   newData,
        Degree: newDegree,
		TwoM: twoM,
    }
}