// metrics_test.go
package main

import "testing"

// ---------- helpers ----------

func almostEqual32(a, b, eps float32) bool {
	if a > b {
		return a-b <= eps
	}
	return b-a <= eps
}

func almostEqual64(a, b, eps float64) bool {
	if a > b {
		return a-b <= eps
	}
	return b-a <= eps
}

// Simple 2-node undirected graph with a single edge (weight 1):
// nodes 0 and 1 connected, symmetric CSR.
//
// N = 2
// Indptr: [0, 1, 2]
// Indices: [1, 0]
// Data: [1, 1]
// Degree: [1, 1]
// TwoM = 2
func buildTwoNodeGraph() *CSR {
	return &CSR{
		N:      2,
		Indptr: []Idx{0, 1, 2},
		Indices: []NodeID{
			1, // neighbor of 0
			0, // neighbor of 1
		},
		Data: []Weight{
			1,
			1,
		},
		Degree: []float32{
			1,
			1,
		},
		TwoM: 2,
	}
}

// ---------- Cosine tests ----------

func TestCosineIdenticalVectors(t *testing.T) {
	a := []float32{1 / 1.4142135, 1 / 1.4142135} // normalized [1,1]
	b := []float32{1 / 1.4142135, 1 / 1.4142135}

	got := Cosine(a, b)
	if !almostEqual32(got, 1.0, 1e-5) {
		t.Fatalf("Cosine(identical) = %v, expected ~1.0", got)
	}
}

func TestCosineOrthogonalVectors(t *testing.T) {
	a := []float32{1, 0}
	b := []float32{0, 1}

	got := Cosine(a, b)
	if !almostEqual32(got, 0.0, 1e-6) {
		t.Fatalf("Cosine(orthogonal) = %v, expected 0.0", got)
	}
}

func TestCosineOppositeVectors(t *testing.T) {
	a := []float32{1, 0}
	b := []float32{-1, 0}

	got := Cosine(a, b)
	if !almostEqual32(got, -1.0, 1e-6) {
		t.Fatalf("Cosine(opposite) = %v, expected -1.0", got)
	}
}

// ---------- deltaQModularity tests ----------

func TestDeltaQModularityZeroAtThreshold(t *testing.T) {
	ki := float32(2.0)
	TotC := float32(4.0)
	TwoM := float32(10.0)

	// For this implementation, ΔQ = 0 when:
	// kin = ki * TotC / TwoM
	threshold := ki * TotC / TwoM

	got := deltaQModularity(threshold, ki, TotC, TwoM)
	if got != 0 {
		t.Fatalf("expected deltaQModularity to be 0 at threshold kin=%v, got %v", threshold, got)
	}
}

func TestDeltaQModularityPositiveWhenKinAboveThreshold(t *testing.T) {
	ki := float32(2.0)
	TotC := float32(4.0)
	TwoM := float32(10.0)

	threshold := ki * TotC / TwoM
	kin := threshold + 0.5

	got := deltaQModularity(kin, ki, TotC, TwoM)
	if got <= 0 {
		t.Fatalf("expected positive deltaQModularity when kin>threshold, got %v", got)
	}
}

func TestDeltaQModularityNegativeWhenKinBelowThreshold(t *testing.T) {
	ki := float32(2.0)
	TotC := float32(4.0)
	TwoM := float32(10.0)

	threshold := ki * TotC / TwoM
	kin := threshold - 0.5

	got := deltaQModularity(kin, ki, TotC, TwoM)
	if got >= 0 {
		t.Fatalf("expected negative deltaQModularity when kin<threshold, got %v", got)
	}
}


// ---------- deltaRBPM tests ----------

func TestDeltaRBPMZeroWhenBalanced(t *testing.T) {
	ki := float32(2.0)
	TotC := float32(4.0)
	TwoM := float32(10.0)
	gamma := float32(1.5)

	kin := gamma * ki * TotC / (2 * TwoM) // note 2*TwoM because formula uses TotC/(TwoM) then /TwoM again
	got := deltaRBPM(kin, ki, TotC, gamma, TwoM)

	// The algebra is (kin - gamma*ki*TotC/TwoM)/TwoM. We want numerator ~0.
	// Here kin = gamma*ki*TotC/(2*TwoM), so it's not exactly 0; but we know
	// sign flips as kin crosses gamma*ki*TotC/TwoM.
	// Better: choose kin = gamma*ki*TotC/TwoM.
	kin = gamma * ki * TotC / TwoM
	got = deltaRBPM(kin, ki, TotC, gamma, TwoM)

	if !almostEqual32(got, 0.0, 1e-6) {
		t.Fatalf("deltaRBPM expected 0, got %v", got)
	}
}

func TestDeltaRBPMSignChange(t *testing.T) {
	ki := float32(2.0)
	TotC := float32(4.0)
	TwoM := float32(10.0)
	gamma := float32(1.5)

	threshold := gamma * ki * TotC / TwoM

	// kin slightly below threshold -> negative
	kinLow := threshold - 0.1
	if deltaRBPM(kinLow, ki, TotC, gamma, TwoM) >= 0 {
		t.Fatalf("expected negative deltaRBPM below threshold, got >= 0")
	}

	// kin slightly above threshold -> positive
	kinHigh := threshold + 0.1
	if deltaRBPM(kinHigh, ki, TotC, gamma, TwoM) <= 0 {
		t.Fatalf("expected positive deltaRBPM above threshold, got <= 0")
	}
}

// ----------
