package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
	"testing"
)

// ============ helpers ============

const leidenFloatEps = 1e-5

func almostEqualFloat32Leiden(a, b float32) bool {
	diff := a - b
	if diff < 0 {
		diff = -diff
	}
	return diff <= leidenFloatEps
}

func sliceEqualInt32Leiden(a, b []int32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func sliceEqualFloat32Leiden(a, b []float32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !almostEqualFloat32Leiden(a[i], b[i]) {
			return false
		}
	}
	return true
}

func splitNonEmptyLeiden(line string) []string {
	return strings.Fields(line)
}

// parse QualityFn from string, e.g. "modularity", "RBPM", "cpm"
func parseQualityFnLeiden(s string) QualityFn {
	s = strings.ToLower(strings.TrimSpace(s))
	switch s {
	case "rbpm":
		return RBPM
	case "cpm":
		return CPM
	case "modularity":
		fallthrough
	default:
		return Modularity
	}
}

// ============ test case structs ============

type LeidenCommunityDetectionTest struct {
	graph      *CSR
	qualityFn  QualityFn
	selfLoop   float32
	gamma      float32
	maxSweeps  int
	maxLevels  int
	capGuess   int
	expPart    Partition
	expStats   []LevelStats
}

type LocalMoveTest struct {
	graph     *CSR
	qualityFn QualityFn
	selfLoop  float32
	gamma     float32
	maxSweeps int
	capGuess  int
	part      Partition
	expMoves  int
	expDeltaQ float32
}

type LocalMoveSweepTest struct {
	graph     *CSR
	qualityFn QualityFn
	selfLoop  float32
	gamma     float32
	capGuess  int
	part      Partition
	order     []int32
	expMoves  int
	expDeltaQ float32
}

type RefinePartitionTest struct {
	graph   *CSR
	part    Partition
	expPart Partition
}

type AggregationTest struct {
	graph   *CSR
	part    Partition
	expCSR  *CSR
}

// ============ top-level tests ============

// ---- LeidenCommunityDetection ----
//This function tests the Leiden Community Detection function
func TestLeidenCommunityDetection(t *testing.T) {
	tests := ReadLeidenCommunityDetectionTests("Tests/LeidenCommunityDetection/")

	for _, test := range tests {
		cfgFn := test.qualityFn
		part, stats := LeidenCommunityDetection(
			test.graph,
			cfgFn,
			test.selfLoop,
			test.gamma,
			test.maxSweeps,
			test.maxLevels,
			test.capGuess,
		)

		if !sliceEqualInt32Leiden(part, test.expPart) {
			t.Fatalf("LeidenCommunityDetection: partition = %v, expected %v",
				part, test.expPart)
		}

		if len(stats) != len(test.expStats) {
			t.Fatalf("LeidenCommunityDetection: len(stats)=%d, expected %d",
				len(stats), len(test.expStats))
		}
		for i := range stats {
			es := test.expStats[i]
			got := stats[i]
			if got.Level != es.Level ||
				got.Moves != es.Moves ||
				got.NumCommunities != es.NumCommunities ||
				!almostEqualFloat32Leiden(float32(got.Quality), float32(es.Quality)) {
				t.Fatalf(
					"LeidenCommunityDetection: LevelStats[%d] = %+v, expected %+v",
					i, got, es,
				)
			}
		}
	}
}

// ---- LocalMove ----
//This function ends up testing the LocalMove function by creating test objects and comparing the actual value to the expected
func TestLocalMove(t *testing.T) {
	tests := ReadLocalMoveTests("Tests/LocalMove/")

	for _, test := range tests {
		cfg := QualityCfg{
			Fn:       test.qualityFn,
			SelfLoop: test.selfLoop,
			Gamma:    test.gamma,
		}
		cs := InitializeCommStatsFromPartition(test.graph, test.part)
		mb := InitializeMoveBuffers(test.capGuess)

		moves, deltaQ := LocalMove(test.graph, cfg, test.part, cs, mb, test.maxSweeps)

		if moves != test.expMoves {
			t.Fatalf("LocalMove: moves=%d, expected %d", moves, test.expMoves)
		}
		if !almostEqualFloat32Leiden(deltaQ, test.expDeltaQ) {
			t.Fatalf("LocalMove: deltaQ=%v, expected %v", deltaQ, test.expDeltaQ)
		}
	}
}

//This function tests the LocalMoveSweep function by comparing the actual value of the test inputs to their expected value
func TestLocalMoveSweepFrom(t *testing.T) {
	tests := ReadLocalMoveSweepTests("Tests/LocalMoveSweep/")
	for _, test := range tests {
		cfg := QualityCfg{
			Fn:       test.qualityFn,
			SelfLoop: test.selfLoop,
			Gamma:    test.gamma,
		}
		cs := InitializeCommStatsFromPartition(test.graph, test.part)
		mb := InitializeMoveBuffers(test.capGuess)

		moves, gains := LocalMoveSweep(test.graph, cfg, test.part, cs, mb, test.order)

		if moves != test.expMoves {
			t.Fatalf("LocalMoveSweep: moves=%d, expected %d", moves, test.expMoves)
		}
		if !almostEqualFloat32Leiden(gains, test.expDeltaQ) {
			t.Fatalf("LocalMoveSweep: gains=%v, expected %v", gains, test.expDeltaQ)
		}
	}
}

// ---- RefinePartition ----
//
// Input (Tests/RefinePartition/input/*.txt):
//   line 1: N
//   line 2: Indptr
//   line 3: Indices
//   line 4: Data              (still required to build CSR; can be all 1s)
//   line 5: Degree (N floats) (can be sums of Data per node; used nowhere here)
//   line 6: TwoM (float)      (unused in RefinePartition)
//   line 7: partition P (N ints)
//
// Output (Tests/RefinePartition/output/*.txt):
//   line 1: refined partition (N ints)
func TestRefinePartition(t *testing.T) {
	tests := ReadRefinePartitionTests("Tests/RefinePartition/")
	for _, test := range tests {
		rb := InitializeRefineBuffers(test.graph, int(test.graph.N))
		ref := RefinePartition(test.graph, test.part, rb)
		if !sliceEqualInt32Leiden(ref, test.expPart) {
			t.Fatalf("RefinePartition: got=%v, expected=%v", ref, test.expPart)
		}
	}
}

// ---- Aggregation ----
func TestAggregation(t *testing.T) {
	tests := ReadAggregationTests("Tests/Aggregation/")

	for _, test := range tests {
		// identity aggregation map to start
		aggMap := make(AggMap, len(test.graph.Degree))
		for i := range aggMap {
			aggMap[i] = int32(i)
		}

		newG := Aggregation(test.graph, test.part, aggMap)
		exp := test.expCSR

		if int(newG.N) != int(exp.N) {
			t.Fatalf("Aggregation: N=%d, expected %d", newG.N, exp.N)
		}
		if len(newG.Indptr) != len(exp.Indptr) ||
			len(newG.Indices) != len(exp.Indices) ||
			len(newG.Data) != len(exp.Data) ||
			len(newG.Degree) != len(exp.Degree) {
			t.Fatalf("Aggregation: CSR length mismatch")
		}

		for i := range newG.Indptr {
			if newG.Indptr[i] != exp.Indptr[i] {
				t.Fatalf("Aggregation: Indptr[%d]=%d, expected %d",
					i, newG.Indptr[i], exp.Indptr[i])
			}
		}
		for i := range newG.Indices {
			if newG.Indices[i] != exp.Indices[i] {
				t.Fatalf("Aggregation: Indices[%d]=%d, expected %d",
					i, newG.Indices[i], exp.Indices[i])
			}
		}
		if !sliceEqualFloat32Leiden(newG.Data, exp.Data) {
			t.Fatalf("Aggregation: Data mismatch")
		}
		if !sliceEqualFloat32Leiden(newG.Degree, exp.Degree) {
			t.Fatalf("Aggregation: Degree mismatch")
		}
		if !almostEqualFloat32Leiden(newG.TwoM, exp.TwoM) {
			t.Fatalf("Aggregation: TwoM=%v, expected %v", newG.TwoM, exp.TwoM)
		}
	}
}

// The following methods end up reading in and creating test objects from various textfiles with different inputs.

func ReadLeidenCommunityDetectionTests(dir string) []LeidenCommunityDetectionTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]LeidenCommunityDetectionTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i].graph,
			tests[i].qualityFn,
			tests[i].selfLoop,
			tests[i].gamma,
			tests[i].maxSweeps,
			tests[i].maxLevels,
			tests[i].capGuess =
			readLeidenInputFile(dir + "input/" + f.Name())
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadLeidenCommunityDetectionTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expPart, tests[i].expStats =
			readLeidenOutputFile(dir + "output/" + f.Name())
	}

	return tests
}

func readLeidenInputFile(file string) (*CSR, QualityFn, float32, float32, int, int, int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// line 1: N
	if !sc.Scan() {
		panic("readLeidenInputFile: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	// Indptr, Indices, Data, Degree, TwoM (CSR)
	indptr := readIdxSliceFromScannerLeiden(sc)
	indices := readNodeIDSliceFromScannerLeiden(sc)
	data := readWeightSliceFromScannerLeiden(sc)
	degree := readFloat32SliceFromScannerLeiden(sc)
	twoM := readSingleFloat32FromScannerLeiden(sc)

	// qualityFn
	if !sc.Scan() {
		panic("readLeidenInputFile: missing qualityFn")
	}
	qFn := parseQualityFnLeiden(sc.Text())

	// selfLoop gamma
	if !sc.Scan() {
		panic("readLeidenInputFile: missing selfLoop/gamma line")
	}
	parts := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	if len(parts) != 2 {
		panic("readLeidenInputFile: selfLoop/gamma line must have 2 floats")
	}
	slf, _ := strconv.ParseFloat(parts[0], 32)
	gam, _ := strconv.ParseFloat(parts[1], 32)

	// maxSweeps maxLevels capGuess
	if !sc.Scan() {
		panic("readLeidenInputFile: missing maxSweeps/maxLevels/capGuess line")
	}
	confParts := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	if len(confParts) != 3 {
		panic("readLeidenInputFile: need 3 ints for maxSweeps/maxLevels/capGuess")
	}
	maxSweeps, _ := strconv.Atoi(confParts[0])
	maxLevels, _ := strconv.Atoi(confParts[1])
	capGuess, _ := strconv.Atoi(confParts[2])

	g := &CSR{
		N:       int32(N),
		Indptr:  indptr,
		Indices: indices,
		Data:    data,
		Degree:  degree,
		TwoM:    twoM,
	}
	return g, qFn, float32(slf), float32(gam), maxSweeps, maxLevels, capGuess
}

func readLeidenOutputFile(file string) (Partition, []LevelStats) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// partition
	if !sc.Scan() {
		panic("readLeidenOutputFile: missing partition line")
	}
	part := readPartitionFromLineLeiden(sc.Text())

	// LevelStats lines
	var stats []LevelStats
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		parts := splitNonEmptyLeiden(line)
		if len(parts) != 4 {
			panic("readLeidenOutputFile: LevelStats line must have 4 fields")
		}
		lv, _ := strconv.Atoi(parts[0])
		qual, _ := strconv.ParseFloat(parts[1], 32)
		nc, _ := strconv.Atoi(parts[2])
		mv, _ := strconv.Atoi(parts[3])

		stats = append(stats, LevelStats{
			Level:         lv,
			Quality:       float64(qual),
			NumCommunities: nc,
			Moves:         mv,
		})
	}
	return part, stats
}

// ---- LocalMove readers ----

func ReadLocalMoveTests(dir string) []LocalMoveTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]LocalMoveTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i] = readLocalMoveInput(dir + "input/" + f.Name())
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadLocalMoveTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		moves, dq := readMovesDeltaQOutput(dir + "output/" + f.Name())
		tests[i].expMoves = moves
		tests[i].expDeltaQ = dq
	}
	return tests
}

func readLocalMoveInput(file string) LocalMoveTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// N
	if !sc.Scan() {
		panic("readLocalMoveInput: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	indptr := readIdxSliceFromScannerLeiden(sc)
	indices := readNodeIDSliceFromScannerLeiden(sc)
	data := readWeightSliceFromScannerLeiden(sc)
	degree := readFloat32SliceFromScannerLeiden(sc)
	twoM := readSingleFloat32FromScannerLeiden(sc)

	// partition P
	if !sc.Scan() {
		panic("readLocalMoveInput: missing partition line")
	}
	part := readPartitionFromLineLeiden(sc.Text())

	// qualityFn
	if !sc.Scan() {
		panic("readLocalMoveInput: missing qualityFn")
	}
	qFn := parseQualityFnLeiden(sc.Text())

	// selfLoop/gamma
	if !sc.Scan() {
		panic("readLocalMoveInput: missing selfLoop/gamma")
	}
	sg := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	if len(sg) != 2 {
		panic("readLocalMoveInput: need 2 floats for selfLoop/gamma")
	}
	slf, _ := strconv.ParseFloat(sg[0], 32)
	gam, _ := strconv.ParseFloat(sg[1], 32)

	// maxSweeps capGuess
	if !sc.Scan() {
		panic("readLocalMoveInput: missing maxSweeps/capGuess")
	}
	mc := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	if len(mc) != 2 {
		panic("readLocalMoveInput: need 2 ints for maxSweeps/capGuess")
	}
	maxSweeps, _ := strconv.Atoi(mc[0])
	capGuess, _ := strconv.Atoi(mc[1])

	return LocalMoveTest{
		graph: &CSR{
			N:       int32(N),
			Indptr:  indptr,
			Indices: indices,
			Data:    data,
			Degree:  degree,
			TwoM:    twoM,
		},
		qualityFn: parseQualityFnLeiden(qFnToStringLeiden(qFn)), // convert to known enum
		selfLoop:  float32(slf),
		gamma:     float32(gam),
		maxSweeps: maxSweeps,
		capGuess:  capGuess,
		part:      part,
	}
}

// helper to serialize QualityFn back to string for parseQualityFnLeiden
func qFnToStringLeiden(q QualityFn) string {
	switch q {
	case RBPM:
		return "rbpm"
	case CPM:
		return "cpm"
	case Modularity:
		fallthrough
	default:
		return "modularity"
	}
}

// ---- LocalMoveSweep readers ----

func ReadLocalMoveSweepTests(dir string) []LocalMoveSweepTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]LocalMoveSweepTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i] = readLocalMoveSweepInput(dir + "input/" + f.Name())
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadLocalMoveSweepTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		mv, dq := readMovesDeltaQOutput(dir + "output/" + f.Name())
		tests[i].expMoves = mv
		tests[i].expDeltaQ = dq
	}
	return tests
}

func readLocalMoveSweepInput(file string) LocalMoveSweepTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	indptr := readIdxSliceFromScannerLeiden(sc)
	indices := readNodeIDSliceFromScannerLeiden(sc)
	data := readWeightSliceFromScannerLeiden(sc)
	degree := readFloat32SliceFromScannerLeiden(sc)
	twoM := readSingleFloat32FromScannerLeiden(sc)

	// partition
	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing partition")
	}
	part := readPartitionFromLineLeiden(sc.Text())

	// order
	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing order")
	}
	order := readInt32SliceFromLineLeiden(sc.Text())

	// qualityFn
	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing qualityFn")
	}
	qFn := parseQualityFnLeiden(sc.Text())

	// selfLoop/gamma
	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing selfLoop/gamma")
	}
	sg := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	if len(sg) != 2 {
		panic("readLocalMoveSweepInput: need 2 floats for selfLoop/gamma")
	}
	slf, _ := strconv.ParseFloat(sg[0], 32)
	gam, _ := strconv.ParseFloat(sg[1], 32)

	// capGuess
	if !sc.Scan() {
		panic("readLocalMoveSweepInput: missing capGuess")
	}
	capGuess, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	return LocalMoveSweepTest{
		graph: &CSR{
			N:       int32(N),
			Indptr:  indptr,
			Indices: indices,
			Data:    data,
			Degree:  degree,
			TwoM:    twoM,
		},
		qualityFn: qFn,
		selfLoop:  float32(slf),
		gamma:     float32(gam),
		capGuess:  capGuess,
		part:      part,
		order:     order,
	}
}

// ---- RefinePartition readers ----

func ReadRefinePartitionTests(dir string) []RefinePartitionTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]RefinePartitionTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i] = readRefinePartitionInput(dir + "input/" + f.Name())
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadRefinePartitionTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expPart = readPartitionFromFileLeiden(dir + "output/" + f.Name())
	}

	return tests
}

func readRefinePartitionInput(file string) RefinePartitionTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readRefinePartitionInput: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	indptr := readIdxSliceFromScannerLeiden(sc)
	indices := readNodeIDSliceFromScannerLeiden(sc)
	data := readWeightSliceFromScannerLeiden(sc)
	degree := readFloat32SliceFromScannerLeiden(sc)
	twoM := readSingleFloat32FromScannerLeiden(sc)

	// partition
	if !sc.Scan() {
		panic("readRefinePartitionInput: missing partition")
	}
	part := readPartitionFromLineLeiden(sc.Text())

	return RefinePartitionTest{
		graph: &CSR{
			N:       int32(N),
			Indptr:  indptr,
			Indices: indices,
			Data:    data,
			Degree:  degree,
			TwoM:    twoM,
		},
		part: part,
	}
}

// ---- Aggregation readers ----

func ReadAggregationTests(dir string) []AggregationTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]AggregationTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i] = readAggregationInput(dir + "input/" + f.Name())
	}

	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadAggregationTests: input/output count mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expCSR = readCSRFromFileLeiden(dir + "output/" + f.Name())
	}
	return tests
}

func readAggregationInput(file string) AggregationTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readAggregationInput: missing N")
	}
	N, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	indptr := readIdxSliceFromScannerLeiden(sc)
	indices := readNodeIDSliceFromScannerLeiden(sc)
	data := readWeightSliceFromScannerLeiden(sc)
	degree := readFloat32SliceFromScannerLeiden(sc)
	twoM := readSingleFloat32FromScannerLeiden(sc)

	if !sc.Scan() {
		panic("readAggregationInput: missing partition")
	}
	part := readPartitionFromLineLeiden(sc.Text())

	return AggregationTest{
		graph: &CSR{
			N:       int32(N),
			Indptr:  indptr,
			Indices: indices,
			Data:    data,
			Degree:  degree,
			TwoM:    twoM,
		},
		part: part,
	}
}

//These are parsers that are supplementary to reading in the different inputs from the files

func readIdxSliceFromScannerLeiden(sc *bufio.Scanner) []Idx {
	if !sc.Scan() {
		panic("readIdxSliceFromScannerLeiden: missing line")
	}
	parts := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	out := make([]Idx, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		out[i] = Idx(v)
	}
	return out
}

func readNodeIDSliceFromScannerLeiden(sc *bufio.Scanner) []NodeID {
	if !sc.Scan() {
		panic("readNodeIDSliceFromScannerLeiden: missing line")
	}
	parts := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	out := make([]NodeID, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		out[i] = NodeID(v)
	}
	return out
}

func readWeightSliceFromScannerLeiden(sc *bufio.Scanner) []Weight {
	if !sc.Scan() {
		panic("readWeightSliceFromScannerLeiden: missing line")
	}
	parts := splitNonEmptyLeiden(strings.TrimSpace(sc.Text()))
	out := make([]Weight, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		out[i] = Weight(v)
	}
	return out
}

func readFloat32SliceFromScannerLeiden(sc *bufio.Scanner) []float32 {
	if !sc.Scan() {
		panic("readFloat32SliceFromScannerLeiden: missing line")
	}
	return readFloat32SliceFromLineLeiden(sc.Text())
}

func readFloat32SliceFromLineLeiden(line string) []float32 {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	out := make([]float32, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		out[i] = float32(v)
	}
	return out
}

func readInt32SliceFromLineLeiden(line string) []int32 {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	out := make([]int32, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 32)
		if err != nil {
			panic(err)
		}
		out[i] = int32(v)
	}
	return out
}

func readSingleFloat32FromScannerLeiden(sc *bufio.Scanner) float32 {
	if !sc.Scan() {
		panic("readSingleFloat32FromScannerLeiden: missing line")
	}
	v, err := strconv.ParseFloat(strings.TrimSpace(sc.Text()), 32)
	if err != nil {
		panic(err)
	}
	return float32(v)
}

func readPartitionFromLineLeiden(line string) Partition {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	P := make(Partition, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 32)
		if err != nil {
			panic(err)
		}
		P[i] = int32(v)
	}
	return P
}

func readPartitionFromFileLeiden(file string) Partition {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	if !sc.Scan() {
		panic("readPartitionFromFileLeiden: empty file")
	}
	return readPartitionFromLineLeiden(sc.Text())
}

func readMovesDeltaQOutput(file string) (int, float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("readMovesDeltaQOutput: missing moves line")
	}
	mv, _ := strconv.Atoi(strings.TrimSpace(sc.Text()))

	if !sc.Scan() {
		panic("readMovesDeltaQOutput: missing deltaQ line")
	}
	dq, err := strconv.ParseFloat(strings.TrimSpace(sc.Text()), 32)
	if err != nil {
		panic(err)
	}
	return mv, float32(dq)
}

func readCSRFromFileLeiden(file string) *CSR {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	// Indptr
	if !sc.Scan() {
		panic("readCSRFromFileLeiden: missing Indptr")
	}
	indptr := readIdxSliceFromLineLeiden(sc.Text())

	// Indices
	if !sc.Scan() {
		panic("readCSRFromFileLeiden: missing Indices")
	}
	indices := readNodeIDSliceFromLineLeiden(sc.Text())

	// Data
	if !sc.Scan() {
		panic("readCSRFromFileLeiden: missing Data")
	}
	data := readWeightSliceFromLineLeiden(sc.Text())

	// Degree
	if !sc.Scan() {
		panic("readCSRFromFileLeiden: missing Degree")
	}
	degree := readFloat32SliceFromLineLeiden(sc.Text())

	// TwoM
	if !sc.Scan() {
		panic("readCSRFromFileLeiden: missing TwoM")
	}
	twoMVal, err := strconv.ParseFloat(strings.TrimSpace(sc.Text()), 32)
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

func readIdxSliceFromLineLeiden(line string) []Idx {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	out := make([]Idx, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		out[i] = Idx(v)
	}
	return out
}

func readNodeIDSliceFromLineLeiden(line string) []NodeID {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	out := make([]NodeID, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseInt(p, 10, 64)
		if err != nil {
			panic(err)
		}
		out[i] = NodeID(v)
	}
	return out
}

func readWeightSliceFromLineLeiden(line string) []Weight {
	parts := splitNonEmptyLeiden(strings.TrimSpace(line))
	out := make([]Weight, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		out[i] = Weight(v)
	}
	return out
}
