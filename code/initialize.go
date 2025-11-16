package main

func InitializeDenseRows() DenseRows {
	// read in csv
	var newData DenseRows

	// newData.N = 
	// newData.D = 
	// newData.Data = 

	return newData
}

func InitializePartition(graph *CSR) Partition {
	par := make(Partition,int(graph.N))
	for i := int32(0); i < graph.N; i++ {
		par[i] = i
	}

	return par
}

// every node in it's own community
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

func InitializeMoveBuffers(capGuess int) *MoveBuffers {
	mb := &MoveBuffers{
		CommIDs: make([]int32,0	,capGuess),
		CommWts: make([]float32,0,capGuess),
		Seen: make(map[int32]int,capGuess),
		Order: make([]int,0,capGuess),
	}
	return mb
}

func (mb *MoveBuffers) Reset() {
    mb.CommIDs = mb.CommIDs[:0]
    mb.CommWts = mb.CommWts[:0]
    mb.Order   = mb.Order[:0]

    for k := range mb.Seen { delete(mb.Seen, k) }
}

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