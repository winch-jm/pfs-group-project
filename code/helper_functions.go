package main
import (
  "gonum.org/v1/gonum/mat"
  "gonum.org/v1/gonum/stat"
)
func PCA(X *mat.Dense, k int) (*mat.Dense, []float64) {
	var pc stat.PC
    ok := pc.PrincipalComponents(X, nil) // centers columns
    if !ok { panic("PCA failed") }

    // Eigenvalues (variance)
    vals := pc.Vars(nil)

    // Project to top-k PCs
    n, _ := X.Dims()
    scores := mat.NewDense(n, k, nil)
    pc.Scores(scores, X, nil)

    return scores, vals[:k]
}

func normalizeExpression(c Cell){
	
}