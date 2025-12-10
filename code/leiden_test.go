// leiden_test.go
package main

import (
	"testing"
)

// ---------- Helpers ----------

// check that i and j are in same community
func sameCommunity(p Partition, i, j int32) bool {
	return p[i] == p[j]
}

// ---------- LeidenCommunityDetection tests ----------

func TestLeidenCommunityDetectionTwoPairs(t *testing.T) {
	g := buildTwoPairGraph()

	cfgFn := Modularity
	selfLoop := float32(0.0)
	gamma := float32(1.0)
	maxSweeps := 10
	maxLevels := 5
	capGuess := 8

	P, stats := LeidenCommunityDetection(g, cfgFn, selfLoop, gamma, maxSweeps, maxLevels, capGuess)

	if len(stats) == 0 {
		t.Fatalf("expected at least one level of stats from LeidenCommunityDetection")
	}

	if CountCommunities(P) != 2 {
		t.Fatalf("expected 2 communities, got %d (P=%v)", CountCommunities(P), P)
	}

	// Expect {0,1} in same community, {2,3} in same community, and the two groups different
	if !sameCommunity(P, 0, 1) {
		t.Fatalf("expected nodes 0 and 1 in same community, got P=%v", P)
	}
	if !sameCommunity(P, 2, 3) {
		t.Fatalf("expected nodes 2 and 3 in same community, got P=%v", P)
	}
	if sameCommunity(P, 0, 2) {
		t.Fatalf("expected (0,1) and (2,3) to be in different communities, got P=%v", P)
	}
}

// ---------- LocalMove tests ----------

func TestLocalMoveMergesCommunities(t *testing.T) {
	g := buildTwoPairGraph()

	P := InitializePartition(g) // initially each node its own community
	cs := InitializeCommStatsFromPartition(g, P)
	mb := InitializeMoveBuffers(8)

	cfg := QualityCfg{
		Fn:       Modularity,
		Gamma:    1.0,
		SelfLoop: 0.0,
	}

	moves, gain := LocalMove(g, cfg, P, cs, mb, 10)

	if moves <= 0 {
		t.Fatalf("expected at least one move in LocalMove, got %d", moves)
	}
	if gain <= 0 {
		t.Fatalf("expected positive total gain in LocalMove, got %f", gain)
	}

	if CountCommunities(P) > 2 {
		t.Fatalf("expected LocalMove to reduce communities to at most 2, got %d (P=%v)",
			CountCommunities(P), P)
	}
}

func TestLocalMoveSweepSinglePass(t *testing.T) {
	g := buildTwoPairGraph()

	P := InitializePartition(g)
	cs := InitializeCommStatsFromPartition(g, P)
	mb := InitializeMoveBuffers(8)

	order := []int32{0, 1, 2, 3}
	cfg := QualityCfg{
		Fn:       Modularity,
		Gamma:    1.0,
		SelfLoop: 0.0,
	}

	initialComms := CountCommunities(P)
	moves, gain := LocalMoveSweep(g, cfg, P, cs, mb, order)

	if moves < 0 {
		t.Fatalf("moves should never be negative, got %d", moves)
	}
	if gain < 0 {
		t.Fatalf("total gain should not be negative, got %f", gain)
	}

	finalComms := CountCommunities(P)
	if finalComms > initialComms {
		t.Fatalf("expected communities to stay same or decrease, initial=%d final=%d",
			initialComms, finalComms)
	}
}

// ---------- RefinePartition tests ----------

func TestRefinePartitionSplitsDisconnectedComponents(t *testing.T) {
	g := buildComponentPlusIsolatedGraph()

	// All nodes initially in the same community label 0
	P := Partition{0, 0, 0, 0}

	rb := InitializeRefineBuffers(g, 8)
	refined := RefinePartition(g, P, rb)

	if CountCommunities(refined) != 3 {
		t.Fatalf("expected 3 refined communities, got %d (refined=%v)",
			CountCommunities(refined), refined)
	}

	// 0 and 1 should still be together
	if !sameCommunity(refined, 0, 1) {
		t.Fatalf("expected 0 and 1 in same refined community, got %v", refined)
	}
	// 2 and 3 isolated from 0/1 and each other (in this graph)
	if sameCommunity(refined, 0, 2) || sameCommunity(refined, 1, 2) {
		t.Fatalf("expected node 2 to be in different community than 0,1; refined=%v", refined)
	}
	if sameCommunity(refined, 0, 3) || sameCommunity(refined, 1, 3) {
		t.Fatalf("expected node 3 to be in different community than 0,1; refined=%v", refined)
	}
}

// ---------- Aggregation tests ----------

func TestAggregationBuildsSupernodeGraph(t *testing.T) {
	g := buildTwoPairGraph()

	// Partition: {0,1} -> comm 0, {2,3} -> comm 1
	part := Partition{0, 0, 1, 1}

	aggMap := make(AggMap, len(part))
	for i := range aggMap {
		aggMap[i] = int32(i)
	}

	newG := Aggregation(g, part, aggMap)

	// Aggregated graph should have 2 supernodes (communities)
	if newG.N != 2 {
		t.Fatalf("expected aggregated graph N=2, got %d", newG.N)
	}

	// Make sure aggregationMap was updated to community ids
	for i, cid := range aggMap {
		if cid != part[i] {
			t.Fatalf("aggregationMap[%d] = %d, expected %d", i, cid, part[i])
		}
	}

	// Expect self-loops on both communities and no inter-community edges
	if len(newG.Indices) != 2 || len(newG.Data) != 2 {
		t.Fatalf("expected 2 edges (two self-loops), got len(Indices)=%d len(Data)=%d",
			len(newG.Indices), len(newG.Data))
	}

	// Check that both edges are self-loops with weight 1.0
	selfLoops := make(map[int32]float32)
	for u := int32(0); u < newG.N; u++ {
		for idx := newG.Indptr[u]; idx < newG.Indptr[u+1]; idx++ {
			v := newG.Indices[idx]
			w := newG.Data[idx]
			if u != v {
				t.Fatalf("expected only self-loops in aggregated graph, found edge %d->%d", u, v)
			}
			selfLoops[u] += w
		}
	}

	if !float32AlmostEqual(selfLoops[0], 1.0, 1e-5) {
		t.Fatalf("expected self-loop weight 1.0 on supernode 0, got %f", selfLoops[0])
	}
	if !float32AlmostEqual(selfLoops[1], 1.0, 1e-5) {
		t.Fatalf("expected self-loop weight 1.0 on supernode 1, got %f", selfLoops[1])
	}

	// Degree should reflect those self-loops
	if len(newG.Degree) != int(newG.N) {
		t.Fatalf("Degree length mismatch: got %d, expected %d", len(newG.Degree), newG.N)
	}
	if !float32AlmostEqual(newG.Degree[0], 1.0, 1e-5) {
		t.Fatalf("expected Degree[0]=1.0, got %f", newG.Degree[0])
	}
	if !float32AlmostEqual(newG.Degree[1], 1.0, 1e-5) {
		t.Fatalf("expected Degree[1]=1.0, got %f", newG.Degree[1])
	}

	// TwoM should be sum of degrees
	expectedTwoM := newG.Degree[0] + newG.Degree[1]
	if !float32AlmostEqual(newG.TwoM, expectedTwoM, 1e-5) {
		t.Fatalf("TwoM mismatch: got %f, expected %f", newG.TwoM, expectedTwoM)
	}
}
