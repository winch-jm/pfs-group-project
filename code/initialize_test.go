// init_test.go
package main

import "testing"

// ---------- helpers ----------

func newCSRWithDegrees(degrees []float32) *CSR {
	return &CSR{
		N:      int32(len(degrees)),
		Degree: degrees,
	}
}

// ---------- InitializePartition tests ----------

func TestInitializePartitionEachNodeOwnCommunity(t *testing.T) {
	g := &CSR{N: 4}
	P := InitializePartition(g)

	if len(P) != int(g.N) {
		t.Fatalf("len(P) = %d, expected %d", len(P), g.N)
	}
	for i := int32(0); i < g.N; i++ {
		if P[i] != i {
			t.Fatalf("P[%d] = %d, expected %d", i, P[i], i)
		}
	}
}

// ---------- InitializeCommStats tests ----------

func TestInitializeCommStatsMatchesDegrees(t *testing.T) {
	deg := []float32{1.0, 2.0, 3.0}
	g := newCSRWithDegrees(deg)

	cs := InitializeCommStats(g)
	if len(cs.Tot) != len(deg) {
		t.Fatalf("len(cs.Tot) = %d, expected %d", len(cs.Tot), len(deg))
	}
	if len(cs.Size) != len(deg) {
		t.Fatalf("len(cs.Size) = %d, expected %d", len(cs.Size), len(deg))
	}
	if len(cs.In) != len(deg) {
		t.Fatalf("len(cs.In) = %d, expected %d", len(cs.In), len(deg))
	}

	for i := range deg {
		if cs.Tot[i] != deg[i] {
			t.Fatalf("cs.Tot[%d] = %v, expected %v", i, cs.Tot[i], deg[i])
		}
		if cs.Size[i] != 1 {
			t.Fatalf("cs.Size[%d] = %d, expected 1", i, cs.Size[i])
		}
		if cs.In[i] != 0 {
			t.Fatalf("cs.In[%d] = %v, expected 0", i, cs.In[i])
		}
	}
}

// ---------- InitializeCommStatsFromPartition tests ----------

func TestInitializeCommStatsFromPartitionAggregatesByCommunity(t *testing.T) {
	// 4 nodes, degrees: [1, 2, 3, 4]
	g := newCSRWithDegrees([]float32{1, 2, 3, 4})

	// Partition: 0,1 -> comm 0; 2,3 -> comm 1
	P := Partition{0, 0, 1, 1}

	cs := InitializeCommStatsFromPartition(g, P)

	if len(cs.Tot) != 2 {
		t.Fatalf("len(cs.Tot) = %d, expected 2", len(cs.Tot))
	}
	if len(cs.Size) != 2 {
		t.Fatalf("len(cs.Size) = %d, expected 2", len(cs.Size))
	}
	if cs.In != nil {
		t.Fatalf("expected cs.In to be nil, got non-nil")
	}

	// Tot[0] = 1+2 = 3, Tot[1] = 3+4 = 7
	if cs.Tot[0] != 3 {
		t.Fatalf("cs.Tot[0] = %v, expected 3", cs.Tot[0])
	}
	if cs.Tot[1] != 7 {
		t.Fatalf("cs.Tot[1] = %v, expected 7", cs.Tot[1])
	}
	if cs.Size[0] != 2 {
		t.Fatalf("cs.Size[0] = %d, expected 2", cs.Size[0])
	}
	if cs.Size[1] != 2 {
		t.Fatalf("cs.Size[1] = %d, expected 2", cs.Size[1])
	}
}

// ---------- CommStats.EnsureCapacity tests ----------

func TestCommStatsEnsureCapacityNoChangeWhenEnough(t *testing.T) {
	cs := &CommStats{
		Tot:  make([]float32, 3),
		In:   make([]float32, 3),
		Size: make([]int32, 3),
	}

	cs.EnsureCapacity(1) // already len=3

	if len(cs.Tot) != 3 || len(cs.Size) != 3 || len(cs.In) != 3 {
		t.Fatalf("EnsureCapacity should not change lengths when capacity sufficient")
	}
}

func TestCommStatsEnsureCapacityGrowsAndPreservesData(t *testing.T) {
	cs := &CommStats{
		Tot:  []float32{1, 2},
		In:   []float32{10, 20},
		Size: []int32{3, 4},
	}

	cs.EnsureCapacity(4) // need length >= 5

	if len(cs.Tot) < 5 {
		t.Fatalf("expected cs.Tot length >= 5, got %d", len(cs.Tot))
	}
	if len(cs.Size) < 5 {
		t.Fatalf("expected cs.Size length >= 5, got %d", len(cs.Size))
	}
	if len(cs.In) < 5 {
		t.Fatalf("expected cs.In length >= 5, got %d", len(cs.In))
	}

	// Original data should be preserved
	if cs.Tot[0] != 1 || cs.Tot[1] != 2 {
		t.Fatalf("cs.Tot not preserved after EnsureCapacity: %v", cs.Tot[:2])
	}
	if cs.In[0] != 10 || cs.In[1] != 20 {
		t.Fatalf("cs.In not preserved after EnsureCapacity: %v", cs.In[:2])
	}
	if cs.Size[0] != 3 || cs.Size[1] != 4 {
		t.Fatalf("cs.Size not preserved after EnsureCapacity: %v", cs.Size[:2])
	}
}

func TestCommStatsEnsureCapacityHandlesNilIn(t *testing.T) {
	cs := &CommStats{
		Tot:  []float32{1},
		In:   nil,
		Size: []int32{1},
	}

	cs.EnsureCapacity(3)

	if len(cs.Tot) < 4 || len(cs.Size) < 4 {
		t.Fatalf("expected Tot/Size length >= 4 after EnsureCapacity")
	}
	// In should still be nil (per code behavior)
	if cs.In != nil {
		t.Fatalf("expected cs.In to remain nil when originally nil")
	}
}

// ---------- InitializeMoveBuffers tests ----------

func TestInitializeMoveBuffersAllocatesCorrectCaps(t *testing.T) {
	capGuess := 10
	mb := InitializeMoveBuffers(capGuess)

	if mb == nil {
		t.Fatalf("InitializeMoveBuffers returned nil")
	}
	if len(mb.CommIDs) != 0 || len(mb.CommWts) != 0 || len(mb.Order) != 0 {
		t.Fatalf("expected empty slices from InitializeMoveBuffers")
	}
	if cap(mb.CommIDs) < capGuess || cap(mb.CommWts) < capGuess || cap(mb.Order) < capGuess {
		t.Fatalf("expected capacities >= capGuess, got CommIDs=%d, CommWts=%d, Order=%d",
			cap(mb.CommIDs), cap(mb.CommWts), cap(mb.Order))
	}
	if len(mb.Seen) != 0 {
		t.Fatalf("expected empty Seen map, got len=%d", len(mb.Seen))
	}
}

// ---------- MoveBuffers.Reset tests ----------

func TestMoveBuffersResetClearsState(t *testing.T) {
	mb := InitializeMoveBuffers(4)

	mb.Add(1, 2.5)
	mb.Add(2, 1.0)

	if len(mb.CommIDs) == 0 || len(mb.CommWts) == 0 || len(mb.Seen) == 0 {
		t.Fatalf("expected non-empty buffers before Reset")
	}

	mb.Reset()

	if len(mb.CommIDs) != 0 || len(mb.CommWts) != 0 || len(mb.Order) != 0 {
		t.Fatalf("expected empty slices after Reset")
	}
	if len(mb.Seen) != 0 {
		t.Fatalf("expected empty Seen map after Reset")
	}
}

// ---------- MoveBuffers.Add tests ----------

func TestMoveBuffersAddNewAndExisting(t *testing.T) {
	mb := InitializeMoveBuffers(4)

	mb.Add(1, 2.0)
	if len(mb.CommIDs) != 1 || len(mb.CommWts) != 1 {
		t.Fatalf("expected one entry after first Add, got CommIDs=%d, CommWts=%d",
			len(mb.CommIDs), len(mb.CommWts))
	}
	if mb.CommIDs[0] != 1 || mb.CommWts[0] != 2.0 {
		t.Fatalf("unexpected first Add state: CommIDs=%v, CommWts=%v", mb.CommIDs, mb.CommWts)
	}

	// Add weight to same community -> no new ID, weight accumulates
	mb.Add(1, 3.0)
	if len(mb.CommIDs) != 1 || len(mb.CommWts) != 1 {
		t.Fatalf("expected still one entry after second Add to same comm, CommIDs=%d, CommWts=%d",
			len(mb.CommIDs), len(mb.CommWts))
	}
	if mb.CommWts[0] != 5.0 {
		t.Fatalf("expected accumulated weight 5.0, got %v", mb.CommWts[0])
	}

	// Add new community
	mb.Add(2, 1.5)
	if len(mb.CommIDs) != 2 || len(mb.CommWts) != 2 {
		t.Fatalf("expected two entries after adding second community, CommIDs=%d, CommWts=%d",
			len(mb.CommIDs), len(mb.CommWts))
	}
}

// ---------- InitializeRefineBuffers tests ----------

func TestInitializeRefineBuffersShape(t *testing.T) {
	g := &CSR{N: 5}
	capGuess := 8

	rb := InitializeRefineBuffers(g, capGuess)

	if rb == nil {
		t.Fatalf("InitializeRefineBuffers returned nil")
	}
	if len(rb.Queue) != 0 {
		t.Fatalf("expected empty Queue on init, got len=%d", len(rb.Queue))
	}
	if cap(rb.Queue) < capGuess {
		t.Fatalf("expected Queue cap >= %d, got %d", capGuess, cap(rb.Queue))
	}
	if len(rb.Visited) != int(g.N) {
		t.Fatalf("expected Visited length %d, got %d", g.N, len(rb.Visited))
	}
	if rb.Stamp != 0 {
		t.Fatalf("expected Stamp=0 at init, got %d", rb.Stamp)
	}
}

// ---------- RefineBuffers.Reset tests ----------

func TestRefineBuffersResetClearsQueueAndIncrementsStamp(t *testing.T) {
	g := &CSR{N: 3}
	rb := InitializeRefineBuffers(g, 4)

	rb.Queue = append(rb.Queue, 0, 1)
	oldStamp := rb.Stamp

	rb.Reset()

	if len(rb.Queue) != 0 {
		t.Fatalf("expected empty Queue after Reset, got len=%d", len(rb.Queue))
	}
	if rb.Stamp != oldStamp+1 {
		t.Fatalf("expected Stamp incremented by 1, old=%d new=%d", oldStamp, rb.Stamp)
	}
	// We intentionally do not clear Visited, so we do NOT assert on it here.
}

// ---------- CountCommunities tests ----------

func TestCountCommunitiesAllDistinct(t *testing.T) {
	P := Partition{0, 1, 2, 3}
	if got := CountCommunities(P); got != 4 {
		t.Fatalf("CountCommunities(%v) = %d, expected 4", P, got)
	}
}

func TestCountCommunitiesWithDuplicates(t *testing.T) {
	P := Partition{0, 1, 0, 2, 2, 1}
	if got := CountCommunities(P); got != 3 {
		t.Fatalf("CountCommunities(%v) = %d, expected 3", P, got)
	}
}
