package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
	"testing"
	"container/heap"
	//"fmt"
)

// ==========================
// Helpers for float equality
// ==========================

const floatEps = 1e-5

func floatsAlmostEqual(a, b float32) bool {
	diff := a - b
	if diff < 0 {
		diff = -diff
	}
	return diff <= floatEps
}

func slicesAlmostEqual(a, b []float32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !floatsAlmostEqual(a[i], b[i]) {
			return false
		}
	}
	return true
}

// ==========================
// Test case structs
// ==========================

// MaxHeap: push a bunch of items, then pop some
type MaxHeapTest struct {
	items          []IndexScore // initial items to push
	numPop         int          // how many pops to perform
	expectedIdx    []int32      // expected indices popped in order
	expectedSim    []float32    // expected similarities popped in order
}

// Similarities: DenseRows -> adjacency matrix (top-k)
type SimilaritiesTest struct {
	N, D, K   int
	data      []float32
	expected  [][]float32
}

// GraphCreation: adjacency matrix -> CSR
type GraphCreationTest struct {
	N        int
	adj      [][]float32
	expected *CSR
}

// WeightedKNN: DenseRows + k -> CSR
type WeightedKNNTest struct {
	N, D, K  int
	data     []float32
	expected *CSR
}

// SimilaritiesParallel: DenseRows + k + numWorkers -> adjacency matrix
type SimilaritiesParallelTest struct {
	N, D, K      int
	numWorkers   int
	data         []float32
	expected     [][]float32
}

// ==========================
// Top-level tests
// ==========================

func TestMaxHeap(t *testing.T) {
	tests := ReadMaxHeapTests("Tests/MaxHeap/")

	for _, test := range tests {
		// Build heap from items
		h := &MaxHeap{}
		for i := range test.items {
			item := &test.items[i]
			heap.Push(h, item)
		}
		var gotIdx []int32
		var gotSim []float32
		for i := 0; i < test.numPop && h.Len() > 0; i++ {
			item := heap.Pop(h).(*IndexScore)
			gotIdx = append(gotIdx, item.index)
			gotSim = append(gotSim, item.similarity)
		}

		if len(gotIdx) != len(test.expectedIdx) ||
			len(gotSim) != len(test.expectedSim) {
			t.Fatalf("MaxHeap: popped lengths mismatch")
		}
		for i := range gotIdx {
			if gotIdx[i] != test.expectedIdx[i] {
				t.Fatalf("MaxHeap: popped index[%d] = %d, expected %d",
					i, gotIdx[i], test.expectedIdx[i])
			}
			if !floatsAlmostEqual(gotSim[i], test.expectedSim[i]) {
				t.Fatalf("MaxHeap: popped sim[%d] = %v, expected %v",
					i, gotSim[i], test.expectedSim[i])
			}
		}
	}
}

func TestSimilarities(t *testing.T) {
	tests := ReadSimilaritiesTests("Tests/Similarities/")

	for _, test := range tests {
		dr := DenseRows{
			N:    test.N,
			D:    test.D,
			Data: test.data,
		}

		adj := Similarities(dr, test.K)
		if len(adj) != len(test.expected) {
			t.Fatalf("Similarities: row count mismatch, got %d, expected %d",
				len(adj), len(test.expected))
		}
		for i := range adj {
			if len(adj[i]) != len(test.expected[i]) {
				t.Fatalf("Similarities: col count mismatch in row %d, got %d, expected %d",
					i, len(adj[i]), len(test.expected[i]))
			}
			for j := range adj[i] {
				if !floatsAlmostEqual(adj[i][j], test.expected[i][j]) {
					t.Fatalf("Similarities: adj[%d][%d] = %v, expected %v",
						i, j, adj[i][j], test.expected[i][j])
				}
			}
		}
	}
}

func TestGraphCreation(t *testing.T) {
	tests := ReadGraphCreationTests("Tests/GraphCreation/")

	for _, test := range tests {
		// DenseRows only used for N in GraphCreation; D/Data are irrelevant.
		dr := DenseRows{N: test.N, D: 0, Data: nil}
		csr := GraphCreation(test.adj, dr)

		exp := test.expected

		if int(csr.N) != test.N {
			t.Fatalf("GraphCreation: N = %d, expected %d", csr.N, test.N)
		}
		if len(csr.Indptr) != len(exp.Indptr) ||
			len(csr.Indices) != len(exp.Indices) ||
			len(csr.Data) != len(exp.Data) ||
			len(csr.Degree) != len(exp.Degree) {
			t.Fatalf("GraphCreation: length mismatch in CSR arrays")
		}

		for i := range csr.Indptr {
			if int(csr.Indptr[i]) != int(exp.Indptr[i]) {
				t.Fatalf("GraphCreation: Indptr[%d] = %d, expected %d",
					i, csr.Indptr[i], exp.Indptr[i])
			}
		}
		for i := range csr.Indices {
			if int(csr.Indices[i]) != int(exp.Indices[i]) {
				t.Fatalf("GraphCreation: Indices[%d] = %d, expected %d",
					i, csr.Indices[i], exp.Indices[i])
			}
		}
		if !slicesAlmostEqual(csr.Data, exp.Data) {
			t.Fatalf("GraphCreation: Data mismatch")
		}
		if !slicesAlmostEqual(csr.Degree, exp.Degree) {
			t.Fatalf("GraphCreation: Degree mismatch")
		}
		if !floatsAlmostEqual(csr.TwoM, exp.TwoM) {
			t.Fatalf("GraphCreation: TwoM = %v, expected %v",
				csr.TwoM, exp.TwoM)
		}
	}
}

func TestWeightedKNN(t *testing.T) {
	tests := ReadWeightedKNNTests("Tests/WeightedKNN/")

	for _, test := range tests {
		dr := DenseRows{
			N:    test.N,
			D:    test.D,
			Data: test.data,
		}

		csr := WeightedKNN(dr, test.K)
		exp := test.expected

		if int(csr.N) != test.N {
			t.Fatalf("WeightedKNN: N = %d, expected %d", csr.N, test.N)
		}
		if len(csr.Indptr) != len(exp.Indptr) ||
			len(csr.Indices) != len(exp.Indices) ||
			len(csr.Data) != len(exp.Data) ||
			len(csr.Degree) != len(exp.Degree) {
			t.Fatalf("WeightedKNN: length mismatch in CSR arrays")
		}

		for i := range csr.Indptr {
			if int(csr.Indptr[i]) != int(exp.Indptr[i]) {
				t.Fatalf("WeightedKNN: Indptr[%d] = %d, expected %d",
					i, csr.Indptr[i], exp.Indptr[i])
			}
		}
		for i := range csr.Indices {
			if int(csr.Indices[i]) != int(exp.Indices[i]) {
				t.Fatalf("WeightedKNN: Indices[%d] = %d, expected %d",
					i, csr.Indices[i], exp.Indices[i])
			}
		}
		if !slicesAlmostEqual(csr.Data, exp.Data) {
			t.Fatalf("WeightedKNN: Data mismatch")
		}
		if !slicesAlmostEqual(csr.Degree, exp.Degree) {
			t.Fatalf("WeightedKNN: Degree mismatch")
		}
		if !floatsAlmostEqual(csr.TwoM, exp.TwoM) {
			t.Fatalf("WeightedKNN: TwoM = %v, expected %v",
				csr.TwoM, exp.TwoM)
		}
	}
}

func TestSimilaritiesParallel(t *testing.T) {
	tests := ReadSimilaritiesParallelTests("Tests/SimilaritiesParallel/")

	for _, test := range tests {
		dr := DenseRows{
			N:    test.N,
			D:    test.D,
			Data: test.data,
		}

		adj := SimilaritiesParallel(dr, test.K, test.numWorkers)

		if len(adj) != len(test.expected) {
			t.Fatalf("SimilaritiesParallel: row count mismatch, got %d, expected %d",
				len(adj), len(test.expected))
		}
		for i := range adj {
			if len(adj[i]) != len(test.expected[i]) {
				t.Fatalf("SimilaritiesParallel: col count mismatch in row %d, got %d, expected %d",
					i, len(adj[i]), len(test.expected[i]))
			}
			for j := range adj[i] {
				if !floatsAlmostEqual(adj[i][j], test.expected[i][j]) {
					t.Fatalf("SimilaritiesParallel: adj[%d][%d] = %v, expected %v",
						i, j, adj[i][j], test.expected[i][j])
				}
			}
		}
	}
}

// ==========================
// Readers for test cases
// ==========================

// Assumes a helper:
// func ReadDirectory(dir string) []os.DirEntry
// already exists in your project (like the N-body harness).

// -------- MaxHeap --------
//
// Input format (Tests/MaxHeap/input/*.txt):
//   line 1: nItems
//   next nItems lines: "<index> <similarity>"
//   next line: nPop
//
// Output format:
//   line 1: indices popped (space-separated int32)
//   line 2: similarities popped (space-separated float32)
func ReadMaxHeapTests(dir string) []MaxHeapTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]MaxHeapTest, len(inputFiles))

	for i, f := range inputFiles {
		items, pops := readMaxHeapInput(dir + "input/" + f.Name())
		tests[i].items = items
		tests[i].numPop = pops
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadMaxHeapTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		idxs, sims := readMaxHeapOutput(dir + "output/" + f.Name())
		tests[i].expectedIdx = idxs
		tests[i].expectedSim = sims
	}
	return tests
}

func readMaxHeapInput(file string) ([]IndexScore, int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// nItems
	if !sc.Scan() {
		panic("readMaxHeapInput: missing nItems")
	}
	nItems, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	items := make([]IndexScore, 0, nItems)

	for i := 0; i < nItems; i++ {
		if !sc.Scan() {
			panic("readMaxHeapInput: missing item line")
		}
		line := strings.TrimSpace(sc.Text())
		parts := strings.Split(line, " ")
		if len(parts) != 2 {
			panic("readMaxHeapInput: item line must be '<index> <similarity>'")
		}
		idx64, err := strconv.ParseInt(parts[0], 10, 32)
		if err != nil {
			panic(err)
		}
		sim64, err := strconv.ParseFloat(parts[1], 32)
		if err != nil {
			panic(err)
		}
		items = append(items, IndexScore{
			index:      int32(idx64),
			similarity: float32(sim64),
			heapIdx:    -1,
		})
	}

	// nPop
	if !sc.Scan() {
		panic("readMaxHeapInput: missing nPop")
	}
	nPop, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	return items, nPop
}

func readMaxHeapOutput(file string) ([]int32, []float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readMaxHeapOutput: missing indices line")
	}
	idxLine := strings.TrimSpace(sc.Text())
	idxParts := splitNonEmpty(idxLine)
	idxs := make([]int32, len(idxParts))
	for i, p := range idxParts {
		v, err := strconv.ParseInt(p, 10, 32)
		if err != nil {
			panic(err)
		}
		idxs[i] = int32(v)
	}

	if !sc.Scan() {
		panic("readMaxHeapOutput: missing similarities line")
	}
	simLine := strings.TrimSpace(sc.Text())
	simParts := splitNonEmpty(simLine)
	sims := make([]float32, len(simParts))
	for i, p := range simParts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		sims[i] = float32(v)
	}
	return idxs, sims
}

// -------- Similarities --------
//
// Input (Tests/Similarities/input/*.txt):
//   line 1: N D K
//   next N lines: D floats each (row of data)
//
// Output (Tests/Similarities/output/*.txt):
//   N lines: N floats each (adjacency matrix)
func ReadSimilaritiesTests(dir string) []SimilaritiesTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]SimilaritiesTest, len(inputFiles))

	for i, f := range inputFiles {
		N, D, K, data := readDenseRowsWithK(dir + "input/" + f.Name())
		tests[i].N = N
		tests[i].D = D
		tests[i].K = K
		tests[i].data = data
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadSimilaritiesTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = readMatrix(dir + "output/" + f.Name())
	}
	return tests
}

// -------- GraphCreation --------
//
// Input (Tests/GraphCreation/input/*.txt):
//   line 1: N
//   next N lines: N floats (adjacency matrix)
//
// Output (Tests/GraphCreation/output/*.txt):
//   line 1: Indptr (space-separated ints)
//   line 2: Indices (space-separated ints)
//   line 3: Data   (space-separated floats)
//   line 4: Degree (space-separated floats)
//   line 5: TwoM   (single float)
func ReadGraphCreationTests(dir string) []GraphCreationTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]GraphCreationTest, len(inputFiles))

	for i, f := range inputFiles {
		N, adj := readGraphCreationInput(dir + "input/" + f.Name())
		tests[i].N = N
		tests[i].adj = adj
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadGraphCreationTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = readCSRFromFile(dir + "output/" + f.Name())
	}
	return tests
}

// -------- WeightedKNN --------
//
// Input (Tests/WeightedKNN/input/*.txt):
//   line 1: N D K
//   next N lines: D floats each
//
// Output (Tests/WeightedKNN/output/*.txt):
//   same 5-line CSR format as GraphCreation
func ReadWeightedKNNTests(dir string) []WeightedKNNTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]WeightedKNNTest, len(inputFiles))

	for i, f := range inputFiles {
		N, D, K, data := readDenseRowsWithK(dir + "input/" + f.Name())
		tests[i].N = N
		tests[i].D = D
		tests[i].K = K
		tests[i].data = data
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadWeightedKNNTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = readCSRFromFile(dir + "output/" + f.Name())
	}
	return tests
}

// -------- SimilaritiesParallel --------
//
// Input (Tests/SimilaritiesParallel/input/*.txt):
//   line 1: N D K numWorkers
//   next N lines: D floats each
//
// Output (Tests/SimilaritiesParallel/output/*.txt):
//   N lines: N floats each (adjacency matrix)
func ReadSimilaritiesParallelTests(dir string) []SimilaritiesParallelTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]SimilaritiesParallelTest, len(inputFiles))

	for i, f := range inputFiles {
		N, D, K, numWorkers, data := readDenseRowsWithKAndWorkers(dir + "input/" + f.Name())
		tests[i].N = N
		tests[i].D = D
		tests[i].K = K
		tests[i].numWorkers = numWorkers
		tests[i].data = data
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadSimilaritiesParallelTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = readMatrix(dir + "output/" + f.Name())
	}
	return tests
}

// ==========================
// Low-level read helpers
// ==========================

func splitNonEmpty(line string) []string {
	if line == "" {
		return []string{}
	}
	parts := strings.Fields(line)
	return parts
}

func readDenseRowsWithK(file string) (int, int, int, []float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// N D K
	if !sc.Scan() {
		panic("readDenseRowsWithK: missing N D K line")
	}
	header := splitNonEmpty(strings.TrimSpace(sc.Text()))
	if len(header) != 3 {
		panic("readDenseRowsWithK: N D K line must have 3 ints")
	}
	N, _ := strconv.Atoi(header[0])
	D, _ := strconv.Atoi(header[1])
	K, _ := strconv.Atoi(header[2])

	data := make([]float32, 0, N*D)
	for i := 0; i < N; i++ {
		if !sc.Scan() {
			panic("readDenseRowsWithK: missing row line")
		}
		line := strings.TrimSpace(sc.Text())
		parts := splitNonEmpty(line)
		if len(parts) != D {
			panic("readDenseRowsWithK: row length != D")
		}
		for _, p := range parts {
			v, err := strconv.ParseFloat(p, 32)
			if err != nil {
				panic(err)
			}
			data = append(data, float32(v))
		}
	}

	return N, D, K, data
}

func readDenseRowsWithKAndWorkers(file string) (int, int, int, int, []float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// N D K numWorkers
	if !sc.Scan() {
		panic("readDenseRowsWithKAndWorkers: missing header line")
	}
	header := splitNonEmpty(strings.TrimSpace(sc.Text()))
	if len(header) != 4 {
		panic("readDenseRowsWithKAndWorkers: header must have 4 ints")
	}
	N, _ := strconv.Atoi(header[0])
	D, _ := strconv.Atoi(header[1])
	K, _ := strconv.Atoi(header[2])
	numWorkers, _ := strconv.Atoi(header[3])

	data := make([]float32, 0, N*D)
	for i := 0; i < N; i++ {
		if !sc.Scan() {
			panic("readDenseRowsWithKAndWorkers: missing row line")
		}
		line := strings.TrimSpace(sc.Text())
		parts := splitNonEmpty(line)
		if len(parts) != D {
			panic("readDenseRowsWithKAndWorkers: row length != D")
		}
		for _, p := range parts {
			v, err := strconv.ParseFloat(p, 32)
			if err != nil {
				panic(err)
			}
			data = append(data, float32(v))
		}
	}
	return N, D, K, numWorkers, data
}

func readMatrix(file string) [][]float32 {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	var rows [][]float32
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		parts := splitNonEmpty(line)
		row := make([]float32, len(parts))
		for i, p := range parts {
			v, err := strconv.ParseFloat(p, 32)
			if err != nil {
				panic(err)
			}
			row[i] = float32(v)
		}
		rows = append(rows, row)
	}
	if err := sc.Err(); err != nil {
		panic(err)
	}
	return rows
}

func readGraphCreationInput(file string) (int, [][]float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readGraphCreationInput: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	rows := make([][]float32, 0, N)
	for i := 0; i < N; i++ {
		if !sc.Scan() {
			panic("readGraphCreationInput: missing row")
		}
		line := strings.TrimSpace(sc.Text())
		parts := splitNonEmpty(line)
		row := make([]float32, len(parts))
		for j, p := range parts {
			v, err := strconv.ParseFloat(p, 32)
			if err != nil {
				panic(err)
			}
			row[j] = float32(v)
		}
		rows = append(rows, row)
	}
	return N, rows
}

func readCSRFromFile(file string) *CSR {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// line 1: Indptr (ints)
	if !sc.Scan() {
		panic("readCSRFromFile: missing Indptr")
	}
	indptrLine := strings.TrimSpace(sc.Text())
	indptrParts := splitNonEmpty(indptrLine)
	indptr := make([]Idx, len(indptrParts))
	for i, p := range indptrParts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		indptr[i] = Idx(v)
	}

	// line 2: Indices (ints)
	if !sc.Scan() {
		panic("readCSRFromFile: missing Indices")
	}
	indicesLine := strings.TrimSpace(sc.Text())
	indicesParts := splitNonEmpty(indicesLine)
	indices := make([]NodeID, len(indicesParts))
	for i, p := range indicesParts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		indices[i] = NodeID(v)
	}

	// line 3: Data (floats)
	if !sc.Scan() {
		panic("readCSRFromFile: missing Data")
	}
	dataLine := strings.TrimSpace(sc.Text())
	dataParts := splitNonEmpty(dataLine)
	data := make([]Weight, len(dataParts))
	for i, p := range dataParts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		data[i] = Weight(v)
	}

	// line 4: Degree (floats)
	if !sc.Scan() {
		panic("readCSRFromFile: missing Degree")
	}
	degLine := strings.TrimSpace(sc.Text())
	degParts := splitNonEmpty(degLine)
	degree := make([]float32, len(degParts))
	for i, p := range degParts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		degree[i] = float32(v)
	}

	// line 5: TwoM (float)
	if !sc.Scan() {
		panic("readCSRFromFile: missing TwoM")
	}
	twoMLine := strings.TrimSpace(sc.Text())
	twoMVal, err := strconv.ParseFloat(twoMLine, 32)
	if err != nil {
		panic(err)
	}

	return &CSR{
		N:       int32(len(degree)),
		Indptr:  indptr,
		Indices: indices,
		Data:    data,
		Degree:  degree,
		TwoM:    float32(twoMVal),
	}
}
