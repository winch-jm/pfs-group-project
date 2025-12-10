package main

import (
	"container/heap"
	"math"
	"testing"
)

// ---------- MaxHeap tests ----------

func TestMaxHeapBasicOrdering(t *testing.T) {
	items := []*IndexScore{
		{similarity: 0.1, index: 0},
		{similarity: 0.5, index: 1},
		{similarity: 0.3, index: 2},
	}

	h := MaxHeap{}
	heap.Init(&h)

	for _, it := range items {
		heap.Push(&h, it)
	}

	// After heap operations, heapIdx should match position
	for i, it := range h {
		if it.heapIdx != i {
			t.Fatalf("heapIdx mismatch after push: item index=%d has heapIdx=%d, expected %d", it.index, it.heapIdx, i)
		}
	}

	// Pop should return in descending similarity order: 0.5, 0.3, 0.1
	expected := []float32{0.5, 0.3, 0.1}
	for i, exp := range expected {
		if h.Len() == 0 {
			t.Fatalf("heap empty when expecting value %v at position %d", exp, i)
		}
		val := heap.Pop(&h).(*IndexScore)
		if math.Abs(float64(val.similarity-exp)) > 1e-6 {
			t.Fatalf("pop %d: got similarity %v, expected %v", i, val.similarity, exp)
		}
		if val.heapIdx != -1 {
			t.Fatalf("popped item heapIdx not -1, got %d", val.heapIdx)
		}
	}
	if h.Len() != 0 {
		t.Fatalf("heap not empty after popping all elements, len=%d", h.Len())
	}
}

func TestMaxHeapSwapUpdatesIndices(t *testing.T) {
	a := &IndexScore{similarity: 0.2, index: 0, heapIdx: 0}
	b := &IndexScore{similarity: 0.8, index: 1, heapIdx: 1}

	h := MaxHeap{a, b}
	h.Swap(0, 1)

	if h[0] != b || h[1] != a {
		t.Fatalf("Swap did not swap elements correctly")
	}
	if h[0].heapIdx != 0 || h[1].heapIdx != 1 {
		t.Fatalf("Swap did not update heapIdx correctly: got (%d, %d), expected (0, 1)",
			h[0].heapIdx, h[1].heapIdx)
	}
}

func TestMaxHeapLenAndLess(t *testing.T) {
	h := MaxHeap{
		{similarity: 0.1, index: 0},
		{similarity: 0.9, index: 1},
	}
	if h.Len() != 2 {
		t.Fatalf("Len() = %d, expected 2", h.Len())
	}
	if !h.Less(1, 0) {
		t.Fatalf("Less() should return true when similarity at i is greater than at j")
	}
	if h.Less(0, 1) {
		t.Fatalf("Less() should return false when similarity at i is less than at j")
	}
}

// ---------- Helpers for DenseRows / float comparison ----------

func float32AlmostEqual(a, b, eps float32) bool {
	return math.Abs(float64(a)-float64(b)) <= float64(eps)
}

// ---------- Similarities tests ----------

func TestSimilaritiesTopK(t *testing.T) {
	// Simple 3-point dataset in 2D:
	// v0 = [1, 0]
	// v1 = [1, 0]  (identical to v0)
	// v2 = [0, 1]  (orthogonal)
	//
	// cos(v0, v1) = 1
	// cos(v0, v2) = 0
	// cos(v1, v2) = 0
	//
	// With k=1, for i=0, the top neighbor should be j=1 with similarity 1.
	dr := DenseRows{
		N: 3,
		D: 2,
		Data: []float32{
			1, 0, // v0
			1, 0, // v1
			0, 1, // v2
		},
	}

	k := 1
	adj := Similarities(dr, k)

	if len(adj) != 3 {
		t.Fatalf("adjacency matrix has wrong number of rows: got %d, expected 3", len(adj))
	}
	for i := range adj {
		if len(adj[i]) != 3 {
			t.Fatalf("adj[%d] has wrong length: got %d, expected 3", i, len(adj[i]))
		}
	}

	// Check that v0 and v1 are connected with similarity ~1
	if !float32AlmostEqual(adj[0][1], 1.0, 1e-6) {
		t.Fatalf("adj[0][1] = %v, expected ~1", adj[0][1])
	}
	if !float32AlmostEqual(adj[1][0], 1.0, 1e-6) {
		t.Fatalf("adj[1][0] = %v, expected ~1", adj[1][0])
	}

	// v0 and v2 should not be in top-1 neighbors for each other, so should be 0
	if !float32AlmostEqual(adj[0][2], 0.0, 1e-6) {
		t.Fatalf("adj[0][2] = %v, expected 0", adj[0][2])
	}
	if !float32AlmostEqual(adj[2][0], 0.0, 1e-6) {
		t.Fatalf("adj[2][0] = %v, expected 0", adj[2][0])
	}
}

// Also confirm k > (n-1) is handled gracefully.
func TestSimilaritiesLargeK(t *testing.T) {
	dr := DenseRows{
		N: 2,
		D: 2,
		Data: []float32{
			1, 0,
			1, 0,
		},
	}

	adj := Similarities(dr, 10) // larger than n-1
	if len(adj) != 2 || len(adj[0]) != 2 || len(adj[1]) != 2 {
		t.Fatalf("adjacency matrix dimensions incorrect for large k")
	}

	// Only one pair of distinct nodes, they should have similarity 1 in both directions
	if !float32AlmostEqual(adj[0][1], 1.0, 1e-6) || !float32AlmostEqual(adj[1][0], 1.0, 1e-6) {
		t.Fatalf("unexpected similarities for large k: adj[0][1]=%v, adj[1][0]=%v",
			adj[0][1], adj[1][0])
	}
}

// ---------- GraphCreation tests ----------

func TestGraphCreationFromAdjacency(t *testing.T) {
	// Use the same adjacency as in TestSimilaritiesTopK for a 3-node graph.
	dr := DenseRows{
		N: 3,
		D: 2,
		Data: []float32{
			1, 0,
			1, 0,
			0, 1,
		},
	}
	adj := Similarities(dr, 1)

	csr := GraphCreation(adj, dr)
	if csr == nil {
		t.Fatalf("GraphCreation returned nil CSR")
	}

	if csr.N != int32(dr.N) {
		t.Fatalf("CSR.N = %d, expected %d", csr.N, dr.N)
	}
	if len(csr.Indptr) != dr.N+1 {
		t.Fatalf("CSR.Indptr length = %d, expected %d", len(csr.Indptr), dr.N+1)
	}

	// Expected: two directed edges: 0->1 and 1->0 (other weights zero)
	// So total edges = 2
	if len(csr.Indices) != 2 {
		t.Fatalf("CSR.Indices length = %d, expected 2", len(csr.Indices))
	}
	if len(csr.Data) != 2 {
		t.Fatalf("CSR.Data length = %d, expected 2", len(csr.Data))
	}

	// Degrees: node 0 -> 1, node 1 -> 1, node 2 -> 0
	if len(csr.Degree) != dr.N {
		t.Fatalf("CSR.Degree length = %d, expected %d", len(csr.Degree), dr.N)
	}
	if !float32AlmostEqual(csr.Degree[0], 1.0, 1e-6) ||
		!float32AlmostEqual(csr.Degree[1], 1.0, 1e-6) ||
		!float32AlmostEqual(csr.Degree[2], 0.0, 1e-6) {
		t.Fatalf("CSR.Degree = %v, expected [1, 1, 0]", csr.Degree)
	}

	// TwoM should be sum of degrees
	if !float32AlmostEqual(csr.TwoM, csr.Degree[0]+csr.Degree[1]+csr.Degree[2], 1e-6) {
		t.Fatalf("CSR.TwoM = %v, expected %v", csr.TwoM,
			csr.Degree[0]+csr.Degree[1]+csr.Degree[2])
	}
}

// ---------- WeightedKNN tests ----------

func TestWeightedKNNConsistency(t *testing.T) {
	// Small 3-point example again
	dr := DenseRows{
		N: 3,
		D: 2,
		Data: []float32{
			1, 0,
			1, 0,
			0, 1,
		},
	}
	k := 1

	// Direct computation
	adj := Similarities(dr, k)
	expectedCSR := GraphCreation(adj, dr)

	// Via WeightedKNN
	gotCSR := WeightedKNN(dr, k)

	if gotCSR == nil {
		t.Fatalf("WeightedKNN returned nil CSR")
	}

	// Check basic fields match
	if gotCSR.N != expectedCSR.N {
		t.Fatalf("WeightedKNN N = %d, expected %d", gotCSR.N, expectedCSR.N)
	}
	if len(gotCSR.Indptr) != len(expectedCSR.Indptr) {
		t.Fatalf("WeightedKNN Indptr length = %d, expected %d",
			len(gotCSR.Indptr), len(expectedCSR.Indptr))
	}
	if len(gotCSR.Indices) != len(expectedCSR.Indices) {
		t.Fatalf("WeightedKNN Indices length = %d, expected %d",
			len(gotCSR.Indices), len(expectedCSR.Indices))
	}
	if len(gotCSR.Data) != len(expectedCSR.Data) {
		t.Fatalf("WeightedKNN Data length = %d, expected %d",
			len(gotCSR.Data), len(expectedCSR.Data))
	}
	if len(gotCSR.Degree) != len(expectedCSR.Degree) {
		t.Fatalf("WeightedKNN Degree length = %d, expected %d",
			len(gotCSR.Degree), len(expectedCSR.Degree))
	}

	// Compare contents (convert indices to int for comparison so we don't care about underlying type)
	for i := range gotCSR.Indptr {
		if int(gotCSR.Indptr[i]) != int(expectedCSR.Indptr[i]) {
			t.Fatalf("Indptr[%d] mismatch: got %d, expected %d",
				i, int(gotCSR.Indptr[i]), int(expectedCSR.Indptr[i]))
		}
	}
	for i := range gotCSR.Indices {
		if int(gotCSR.Indices[i]) != int(expectedCSR.Indices[i]) {
			t.Fatalf("Indices[%d] mismatch: got %d, expected %d",
				i, int(gotCSR.Indices[i]), int(expectedCSR.Indices[i]))
		}
	}
	for i := range gotCSR.Data {
		if !float32AlmostEqual(gotCSR.Data[i], expectedCSR.Data[i], 1e-6) {
			t.Fatalf("Data[%d] mismatch: got %v, expected %v",
				i, gotCSR.Data[i], expectedCSR.Data[i])
		}
	}
	for i := range gotCSR.Degree {
		if !float32AlmostEqual(gotCSR.Degree[i], expectedCSR.Degree[i], 1e-6) {
			t.Fatalf("Degree[%d] mismatch: got %v, expected %v",
				i, gotCSR.Degree[i], expectedCSR.Degree[i])
		}
	}

	if !float32AlmostEqual(gotCSR.TwoM, expectedCSR.TwoM, 1e-6) {
		t.Fatalf("TwoM mismatch: got %v, expected %v", gotCSR.TwoM, expectedCSR.TwoM)
	}
}