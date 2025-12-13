// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

// Initialize every node to be in it's own community
func InitializePartition(graph *CSR) Partition {
	par := make(Partition,int(graph.N))
	for i := int32(0); i < graph.N; i++ {
		par[i] = i
	}

	return par
}

// Initially all communities have: 
// - no internal edges, 
// - size 1 
// - Tot (sum of degrees) == degree of single vertex
func InitializeCommStats(graph *CSR) *CommStats {
	cs := &CommStats{
		Tot: make([]float32,graph.N),
		In: make([]float32,graph.N),
		Size: make([]int32,graph.N),
	}

	for i := int32(0); i < graph.N; i++ {
		cs.Tot[i] = graph.Degree[i]
		cs.Size[i] = 1
	}
	return cs
}


// Derive Community Stats from a given Partition (list of community labels)
func InitializeCommStatsFromPartition(g *CSR, P Partition) *CommStats{
	var maxC int32
	for _, c := range P {
		if c > maxC {
			maxC = c
		}
	}

	C := int(maxC) + 1

	cs := &CommStats{
		Tot:make([]float32, C),
		Size:make([]int32, C),
		In:nil,
	}

	for i := int32(0); i < g.N; i++ {
		c := P[i]
		cs.Tot[c] += g.Degree[i]
		cs.Size[c]++
	}

	return cs
}

// Check if CommStats arrays have enough allocated memory:
// if yes, do nothing
// if no, increase memory allocation
func (cs *CommStats) EnsureCapacity(c int32) {
	if int(c) < len(cs.Tot) {
		return // already large enough
	}

	newLen := int(c) + 1

	// Grow Tot
	oldTot := cs.Tot
	cs.Tot = make([]float32, newLen)
	copy(cs.Tot, oldTot)

	// Grow Size
    oldSize := cs.Size
    cs.Size = make([]int32, newLen)
    copy(cs.Size, oldSize)

	// Grow In (optional)
    if cs.In != nil {
        oldIn := cs.In
        cs.In = make([]float32, newLen)
        copy(cs.In, oldIn)
    }


}

// Allocate buffers for local moving (of prescribed size)
func InitializeMoveBuffers(capGuess int) *MoveBuffers {
	mb := &MoveBuffers{
		CommIDs: make([]int32,0	,capGuess),
		CommWts: make([]float32,0,capGuess),
		Seen: make(map[int32]int,capGuess),
		Order: make([]int,0,capGuess),
	}
	return mb
}

// clear out move buffers to reuse
func (mb *MoveBuffers) Reset() {
    mb.CommIDs = mb.CommIDs[:0]
    mb.CommWts = mb.CommWts[:0]
    mb.Order   = mb.Order[:0]

    for k := range mb.Seen { delete(mb.Seen, k) }
}

// Add 
func (mb *MoveBuffers) Add(commID int32, w float32) {
	if idx, ok := mb.Seen[commID]; ok {
		// community already seen --> accumulate weight
		mb.CommWts[idx] += w
		return
	}
	// new community --> append
	idx := len(mb.CommIDs)
	mb.Seen[commID] = idx
	mb.CommIDs = append(mb.CommIDs,commID)
	mb.CommWts = append(mb.CommWts,w)
}

func InitializeRefineBuffers(graph *CSR, capGuess int) *RefineBuffers {
	rb := &RefineBuffers{
		Queue: make([]int32,0,capGuess),
		Visited: make([]int32,graph.N),
		Stamp: 0,
	}
	return rb

}

func (rb *RefineBuffers) Reset() {
	rb.Queue = rb.Queue[:0]
	rb.Stamp++ // new "epoch"
	// Visited is left as-is; stamp logic makes old marks obsolete
}


func CountCommunities(P Partition) int {
    seen := make(map[int32]struct{})
    for _, c := range P {
        seen[c] = struct{}{}
    }
    return len(seen)
}