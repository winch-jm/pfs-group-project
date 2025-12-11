// Authors: Jeff Winchell, Ajay Prabhakar
// Date: 12/09/2025
package main

import (
	// "math"
)

// cosine simlarity between two normalized vectors
func Cosine(a, b []float32) float32 {

	var s float32
	// dot product
	for i := range a {
		s += a[i]*b[i]
	}

	return s
}

// Leiden Metrics

// Not completely accurate modularity calculations, but works for choosing best move
func deltaQModularity(kin, ki, TotC, TwoM float32) float32 {
	// ΔQ = ( k_i,in(c) - k_i * Tot_c / (2m) ) / (2m)
	return (kin - ki*TotC/(TwoM))/TwoM
}

// for sanity checking modularity
func ComputeModularity(g *CSR, P Partition) float64 {

    twoM := float64(g.TwoM)
    if twoM == 0 {
        return 0
    }

    // find max community id
    maxC := int32(0)
    for _, c := range P {
        if c > maxC {
            maxC = c
        }
    }
    C := int(maxC) + 1

    Tot := make([]float64, C)
    In  := make([]float64, C)

    // Tot[c] = sum of degrees in community c
    for i := int32(0); i < g.N; i++ {
        c := P[i]
        Tot[c] += float64(g.Degree[i])
    }

    // In[c] = sum of internal weights in community c
    // We’ll count each undirected edge once by only taking i<j
    for i := int32(0); i < g.N; i++ {
        ci := P[i]
        rowStart := g.Indptr[i]
        rowEnd   := g.Indptr[i+1]
        for idx := rowStart; idx < rowEnd; idx++ {
            j := g.Indices[idx]
            if j < i {
                continue // avoid double-counting
            }
            if P[j] != ci {
                continue
            }
            w := g.Data[idx]
            In[ci] += float64(w)
        }
    }

    var Q float64
    for c := 0; c < C; c++ {
        if Tot[c] == 0 {
            continue
        }
        Q += In[c]/twoM - (Tot[c]*Tot[c])/(4*twoM*twoM)
    }

    return Q
}

func deltaRBPM(kin,  ki, TotC, gamma, TwoM float32) float32 {
	return (kin - gamma*ki*TotC/(TwoM))/TwoM
}

func deltaCPM(kin, gamma float32, cSize int32) float32 {
	return kin - gamma*float32(cSize)
}



