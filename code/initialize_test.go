package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
	"testing"
	"io/fs"
	"path/filepath"
	"fmt"
)

// ========================
// Test struct definitions
// ========================

type InitializePartitionTest struct {
	N       int32
	PExpect Partition
}

type InitCommStatsTest struct {
	graphN   int32
	degrees  []float32
	expected *CommStats
}

type InitCommStatsFromPartTest struct {
	graphN   int32
	degrees  []float32
	partition Partition
	expected *CommStats
}

type EnsureCapacityTest struct {
	initial *CommStats
	c       int32
	expected *CommStats
}

type InitMoveBuffersTest struct {
	capGuess int
	expectedCap int
}

type MoveBuffersAddTest struct {
	capGuess int
	ops      []MoveBuffersOp
	expectedCommIDs []int32
	expectedCommWts []float32
}

type MoveBuffersOp struct {
	CommID int32
	W      float32
}

type MoveBuffersResetTest struct {
	capGuess int
	opsBeforeReset []MoveBuffersOp
	// After reset, all slices empty
}

type InitRefineBuffersTest struct {
	graphN   int32
	capGuess int
	expectedQueueCap int
}

type RefineBuffersResetTest struct {
	graphN   int32
	capGuess int
	initialQueue []int32
}

type CountCommunitiesTest struct {
	P       Partition
	result  int
}

// ========================
// Top-level Test functions
// ========================

// ---------- InitializePartition ----------

func TestInitializePartition(t *testing.T) {
	tests := ReadInitializePartitionTests("Tests/InitializePartition/")

	for _, test := range tests {
		g := &CSR{N: test.N}
		P := InitializePartition(g)

		if len(P) != int(test.N) {
			t.Fatalf("len(P) = %d, expected %d", len(P), test.N)
		}
		for i := int32(0); i < test.N; i++ {
			if P[i] != test.PExpect[i] {
				t.Fatalf("P[%d] = %d, expected %d", i, P[i], test.PExpect[i])
			}
		}
	}
}
// ---------- InitializeCommStats ----------

func TestInitializeCommStats(t *testing.T) {
	tests := ReadInitCommStatsTests("Tests/InitializeCommStats/")
	for _, test := range tests {
		g := &CSR{
			N:      test.graphN,
			Degree: test.degrees,
		}
		cs := InitializeCommStats(g)

		if len(cs.Tot) != len(test.expected.Tot) ||
			len(cs.Size) != len(test.expected.Size) ||
			len(cs.In)  != len(test.expected.In) {
			t.Fatalf("InitializeCommStats: length mismatch")
		}

		for i := range cs.Tot {
			if cs.Tot[i] != test.expected.Tot[i] {
				t.Fatalf("InitializeCommStats: Tot[%d] = %v, expected %v",
					i, cs.Tot[i], test.expected.Tot[i])
			}
			if cs.Size[i] != test.expected.Size[i] {
				t.Fatalf("InitializeCommStats: Size[%d] = %d, expected %d",
					i, cs.Size[i], test.expected.Size[i])
			}
			if cs.In[i] != test.expected.In[i] {
				t.Fatalf("InitializeCommStats: In[%d] = %v, expected %v",
					i, cs.In[i], test.expected.In[i])
			}
		}
	}
}

// ---------- InitializeCommStatsFromPartition ----------

func TestInitializeCommStatsFromPartition(t *testing.T) {
	tests := ReadInitCommStatsFromPartTests("Tests/InitializeCommStatsFromPartition/")
	for _, test := range tests {
		g := &CSR{
			N:      test.graphN,
			Degree: test.degrees,
		}
		cs := InitializeCommStatsFromPartition(g, test.partition)

		if len(cs.Tot) != len(test.expected.Tot) ||
			len(cs.Size) != len(test.expected.Size) {
			t.Fatalf("InitializeCommStatsFromPartition: length mismatch")
		}
		for i := range cs.Tot {
			if cs.Tot[i] != test.expected.Tot[i] {
				t.Fatalf("InitializeCommStatsFromPartition: Tot[%d] = %v, expected %v",
					i, cs.Tot[i], test.expected.Tot[i])
			}
			if cs.Size[i] != test.expected.Size[i] {
				t.Fatalf("InitializeCommStatsFromPartition: Size[%d] = %d, expected %d",
					i, cs.Size[i], test.expected.Size[i])
			}
		}
	}
}

// ---------- CommStats.EnsureCapacity ----------

func TestEnsureCapacity(t *testing.T) {
	tests := ReadEnsureCapacityTests("Tests/EnsureCapacity/")
	for _, test := range tests {
		cs := &CommStats{
			Tot:  append([]float32(nil), test.initial.Tot...),
			In:   nil,
			Size: append([]int32(nil), test.initial.Size...),
		}
		if test.initial.In != nil {
			cs.In = append([]float32(nil), test.initial.In...)
		}

		cs.EnsureCapacity(test.c)

		// compare with expected
		if len(cs.Tot) != len(test.expected.Tot) ||
			len(cs.Size) != len(test.expected.Size) ||
			((cs.In == nil) != (test.expected.In == nil)) ||
			(cs.In != nil && len(cs.In) != len(test.expected.In)) {

			fmt.Println(cs.In, test.expected.In)
			t.Fatalf("EnsureCapacity: length or nil mismatch after EnsureCapacity(c=%d)", test.c)
		}

		for i := range test.expected.Tot {
			if cs.Tot[i] != test.expected.Tot[i] {
				t.Fatalf("EnsureCapacity: Tot[%d] = %v, expected %v", i, cs.Tot[i], test.expected.Tot[i])
			}
		}
		for i := range test.expected.Size {
			if cs.Size[i] != test.expected.Size[i] {
				t.Fatalf("EnsureCapacity: Size[%d] = %d, expected %d", i, cs.Size[i], test.expected.Size[i])
			}
		}
		if cs.In != nil {
			for i := range test.expected.In {
				if cs.In[i] != test.expected.In[i] {
					t.Fatalf("EnsureCapacity: In[%d] = %v, expected %v", i, cs.In[i], test.expected.In[i])
				}
			}
		}
	}
}

// ---------- InitializeMoveBuffers ----------

func TestInitializeMoveBuffers(t *testing.T) {
	tests := ReadInitMoveBuffersTests("Tests/InitializeMoveBuffers/")
	for _, test := range tests {
		mb := InitializeMoveBuffers(test.capGuess)
		if cap(mb.CommIDs) < test.expectedCap ||
			cap(mb.CommWts) < test.expectedCap ||
			cap(mb.Order)  < test.expectedCap {
			t.Fatalf("InitializeMoveBuffers(capGuess=%d): capacities too small", test.capGuess)
		}
		if len(mb.CommIDs) != 0 || len(mb.CommWts) != 0 || len(mb.Order) != 0 || len(mb.Seen) != 0 {
			t.Fatalf("InitializeMoveBuffers: expected empty slices and map")
		}
	}
}

// ---------- MoveBuffers.Add ----------

func TestMoveBuffersAdd(t *testing.T) {
	tests := ReadMoveBuffersAddTests("Tests/MoveBuffersAdd/")
	for _, test := range tests {
		mb := InitializeMoveBuffers(test.capGuess)
		for _, op := range test.ops {
			mb.Add(op.CommID, op.W)
		}
		if len(mb.CommIDs) != len(test.expectedCommIDs) ||
			len(mb.CommWts) != len(test.expectedCommWts) {
			t.Fatalf("MoveBuffers.Add: length mismatch")
		}
		for i := range mb.CommIDs {
			if mb.CommIDs[i] != test.expectedCommIDs[i] {
				t.Fatalf("MoveBuffers.Add: CommIDs[%d] = %d, expected %d",
					i, mb.CommIDs[i], test.expectedCommIDs[i])
			}
			if mb.CommWts[i] != test.expectedCommWts[i] {
				t.Fatalf("MoveBuffers.Add: CommWts[%d] = %v, expected %v",
					i, mb.CommWts[i], test.expectedCommWts[i])
			}
		}
	}
}

// ---------- MoveBuffers.Reset ----------

func TestMoveBuffersReset(t *testing.T) {
	tests := ReadMoveBuffersResetTests("Tests/MoveBuffersReset/")
	for _, test := range tests {
		mb := InitializeMoveBuffers(test.capGuess)
		for _, op := range test.opsBeforeReset {
			mb.Add(op.CommID, op.W)
		}
		mb.Reset()
		if len(mb.CommIDs) != 0 || len(mb.CommWts) != 0 || len(mb.Order) != 0 {
			t.Fatalf("MoveBuffers.Reset: expected empty slices after reset")
		}
		if len(mb.Seen) != 0 {
			t.Fatalf("MoveBuffers.Reset: expected empty Seen map after reset")
		}
	}
}

// ---------- InitializeRefineBuffers ----------

func TestInitializeRefineBuffers(t *testing.T) {
	tests := ReadInitRefineBuffersTests("Tests/InitializeRefineBuffers/")
	for _, test := range tests {
		g := &CSR{N: test.graphN}
		rb := InitializeRefineBuffers(g, test.capGuess)
		if len(rb.Queue) != 0 || rb.Stamp != 0 {
			t.Fatalf("InitializeRefineBuffers: expected empty queue and Stamp=0")
		}
		if cap(rb.Queue) < test.expectedQueueCap {
			t.Fatalf("InitializeRefineBuffers: Queue cap=%d, expected >=%d", cap(rb.Queue), test.expectedQueueCap)
		}
		if len(rb.Visited) != int(g.N) {
			t.Fatalf("InitializeRefineBuffers: Visited length=%d, expected %d", len(rb.Visited), g.N)
		}
	}
}


// ---------- CountCommunities ----------

func TestCountCommunities(t *testing.T) {
	tests := ReadCountCommunitiesTests("Tests/CountCommunities/")
	for _, test := range tests {
		got := CountCommunities(test.P)
		if got != test.result {
			t.Fatalf("CountCommunities(%v) = %d, expected %d", test.P, got, test.result)
		}
	}
}

// ========================
// Readers for test cases
// ========================

// --- InitializePartition ---

// Format input (per file):
// line 1: N
//
// Format output:
// line 1: space-separated partition labels (int32).

func ReadInitializePartitionTests(directory string) []InitializePartitionTest {
	// INPUT FILES
	rawInputs := ReadDirectory(directory + "/input")

	var inputFiles []os.DirEntry
	for _, f := range rawInputs {
		name := f.Name()
		// skip hidden / metadata files (macOS, etc.)
		if strings.HasPrefix(name, ".") || strings.HasPrefix(name, "._") {
			continue
		}
		// only accept .txt
		if filepath.Ext(name) != ".txt" {
			continue
		}
		inputFiles = append(inputFiles, f)
	}

	numFiles := len(inputFiles)
	tests := make([]InitializePartitionTest, numFiles)

	// read N from input files
	for i, inputFile := range inputFiles {
		tests[i].N = ReadInitializePartitionInput(
			directory + "input/" + inputFile.Name(),
		)
	}

	// OUTPUT FILES
	rawOutputs := ReadDirectory(directory + "/output")

	var outputFiles []os.DirEntry
	for _, f := range rawOutputs {
		name := f.Name()
		if strings.HasPrefix(name, ".") || strings.HasPrefix(name, "._") {
			continue
		}
		if filepath.Ext(name) != ".txt" {
			continue
		}
		outputFiles = append(outputFiles, f)
	}

	if len(outputFiles) != numFiles {
		panic("ReadInitializePartitionTests: number of input and output files do not match")
	}

	// read expected partitions
	for i, outputFile := range outputFiles {
		tests[i].PExpect = ReadPartitionFromFile(
			directory + "output/" + outputFile.Name(),
		)
	}

	return tests
}


// --- InitializeCommStats ---

// Input:
//   line 1: N
//   line 2: degrees (space-separated floats)
// Output:
//   line 1: Tot (space-separated floats)
//   line 2: Size (space-separated int32)
//   line 3: In   (space-separated floats)
func ReadInitCommStatsTests(dir string) []InitCommStatsTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]InitCommStatsTest, len(inputFiles))

	for i, f := range inputFiles {
		N, degrees := ReadGraphDegreesInput(dir + "input/" + f.Name())
		tests[i].graphN = N
		tests[i].degrees = degrees
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadInitCommStatsTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		cs := ReadCommStatsFromFile(dir + "output/" + f.Name())
		tests[i].expected = cs
	}
	return tests
}
func ReadInitializePartitionInput(file string) int32 {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	if !sc.Scan() {
		panic("ReadInitializePartitionInput: empty input file")
	}
	line := strings.TrimSpace(sc.Text())
	if line == "" {
		panic("ReadInitializePartitionInput: missing N")
	}

	n64, err := strconv.ParseInt(line, 10, 32)
	if err != nil {
		panic(err)
	}
	return int32(n64)
}
// --- InitializeCommStatsFromPartition ---

// Input:
//   line 1: N
//   line 2: degrees (floats)
//   line 3: partition labels (int32)
// Output: same 3-line format as CommStats above.
func ReadInitCommStatsFromPartTests(dir string) []InitCommStatsFromPartTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]InitCommStatsFromPartTest, len(inputFiles))

	for i, f := range inputFiles {
		N, deg, P := ReadGraphDegreesAndPartitionInput(dir + "input/" + f.Name())
		tests[i].graphN = N
		tests[i].degrees = deg
		tests[i].partition = P
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadInitCommStatsFromPartTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = ReadCommStatsFromFile(dir + "output/" + f.Name())
	}
	return tests
}

// --- EnsureCapacity ---

// Input format:
//   line 1: Tot (floats)
//   line 2: Size (int32)
//   line 3: In or "nil" (floats or literal "nil")
//   line 4: c (int32)
//
// Output format (final CommStats):
//   line 1: Tot
//   line 2: Size
//   line 3: In or "nil"
func ReadEnsureCapacityTests(dir string) []EnsureCapacityTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]EnsureCapacityTest, len(inputFiles))

	for i, f := range inputFiles {
		cs, c := ReadEnsureCapacityInput(dir + "input/" + f.Name())
		tests[i].initial = cs
		tests[i].c = c
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadEnsureCapacityTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		tests[i].expected = ReadCommStatsPartialFromFile(dir + "output/" + f.Name())
	}
	return tests
}

// --- InitializeMoveBuffers ---

// Input:
//   line 1: capGuess
// Output:
//   line 1: expectedCap (we check cap >= expectedCap)
func ReadInitMoveBuffersTests(dir string) []InitMoveBuffersTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]InitMoveBuffersTest, len(inputFiles))

	for i, f := range inputFiles {
		capGuess, err := readIntFromFile(dir + "input/" + f.Name())
		if err != nil {
			panic(err)
		}
		tests[i].capGuess = capGuess
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadInitMoveBuffersTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		expCap, err := readIntFromFile(dir + "output/" + f.Name())
		if err != nil {
			panic(err)
		}
		tests[i].expectedCap = expCap
	}
	return tests
}

// --- MoveBuffers.Add ---

// Input:
//   line 1: capGuess
//   line 2: number of Add operations (k)
//   next k lines: "commID w" (int32 float32)
// Output:
//   line 1: final CommIDs
//   line 2: final CommWts
func ReadMoveBuffersAddTests(dir string) []MoveBuffersAddTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]MoveBuffersAddTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i].capGuess, tests[i].ops = ReadMoveBuffersAddInput(dir + "input/" + f.Name())
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadMoveBuffersAddTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		ids, wts := ReadInt32SliceFromFile(dir + "output/" + f.Name())
		tests[i].expectedCommIDs = ids
		tests[i].expectedCommWts = wts
	}
	return tests
}

// --- MoveBuffers.Reset ---

// Input:
//   line 1: capGuess
//   line 2: k (number of Adds before reset)
//   next k lines: "commID w"
// Output: (just a file to match count; contents unused)
func ReadMoveBuffersResetTests(dir string) []MoveBuffersResetTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]MoveBuffersResetTest, len(inputFiles))

	for i, f := range inputFiles {
		tests[i].capGuess, tests[i].opsBeforeReset = ReadMoveBuffersAddInput(dir + "input/" + f.Name())
	}
	// Outputs exist only to satisfy the harness pattern; nothing is actually read.
	return tests
}

// --- InitializeRefineBuffers ---

// Input:
//   line 1: N
//   line 2: capGuess
// Output:
//   line 1: expectedQueueCap
func ReadInitRefineBuffersTests(dir string) []InitRefineBuffersTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]InitRefineBuffersTest, len(inputFiles))

	for i, f := range inputFiles {
		N, capGuess := ReadRefineBuffersInitInput(dir + "input/" + f.Name())
		tests[i].graphN = N
		tests[i].capGuess = capGuess
	}
	outputFiles := ReadDirectory(dir + "/output")
	if len(outputFiles) != len(inputFiles) {
		panic("ReadInitRefineBuffersTests: input/output mismatch")
	}
	for i, f := range outputFiles {
		expCap, err := readIntFromFile(dir + "output/" + f.Name())
		if err != nil {
			panic(err)
		}
		tests[i].expectedQueueCap = expCap
	}
	return tests
}

// --- RefineBuffers.Reset ---

// Input:
//   line 1: N
//   line 2: capGuess
//   line 3: space-separated queue contents (int32)
// Output: dummy
func ReadRefineBuffersResetTests(dir string) []RefineBuffersResetTest {
	inputFiles := ReadDirectory(dir + "/input")
	tests := make([]RefineBuffersResetTest, len(inputFiles))

	for i, f := range inputFiles {
		N, capGuess, queue := ReadRefineBuffersResetInput(dir + "input/" + f.Name())
		tests[i].graphN = N
		tests[i].capGuess = capGuess
		tests[i].initialQueue = queue
	}
	return tests
}

// --- CountCommunities ---

// Input:
//   line 1: partition labels
// Output:
//   line 1: integer result
func ReadCountCommunitiesTests(dir string) []CountCommunitiesTest {
	// --- INPUT FILES ---

	rawInputs := ReadDirectory(dir + "/input")

	var inputFiles []os.DirEntry
	for _, f := range rawInputs {
		name := f.Name()
		// skip hidden / metadata files (macOS, etc.)
		if strings.HasPrefix(name, ".") || strings.HasPrefix(name, "._") {
			continue
		}
		// only accept .txt files
		if !strings.HasSuffix(name, ".txt") {
			continue
		}
		inputFiles = append(inputFiles, f)
	}

	tests := make([]CountCommunitiesTest, len(inputFiles))

	// read partitions from input files
	for i, f := range inputFiles {
		tests[i].P = ReadPartitionFromFile(dir + "input/" + f.Name())
	}

	// --- OUTPUT FILES ---

	rawOutputs := ReadDirectory(dir + "/output")

	var outputFiles []os.DirEntry
	for _, f := range rawOutputs {
		name := f.Name()
		if strings.HasPrefix(name, ".") || strings.HasPrefix(name, "._") {
			continue
		}
		if !strings.HasSuffix(name, ".txt") {
			continue
		}
		outputFiles = append(outputFiles, f)
	}

	if len(outputFiles) != len(inputFiles) {
		panic("ReadCountCommunitiesTests: input/output mismatch")
	}

	for i, f := range outputFiles {
		val32, err := readInt32FromFile(dir + "output/" + f.Name())
		if err != nil {
			panic(err)
		}
		// assuming CountCommunitiesTest.result is an int
		tests[i].result = int(val32)
	}

	return tests
}
// ========================
// Low-level input helpers
// ========================

func readInt32FromFile(file string) (int32, error) {
	f, err := os.Open(file)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	if !sc.Scan() {
		return 0, sc.Err()
	}
	line := strings.TrimSpace(sc.Text())
	n64, err := strconv.ParseInt(line, 10, 32)
	return int32(n64), err
}

func readIntFromFile(file string) (int, error) {
	f, err := os.Open(file)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	if !sc.Scan() {
		return 0, sc.Err()
	}
	line := strings.TrimSpace(sc.Text())
	n64, err := strconv.ParseInt(line, 10, 64)
	return int(n64), err
}

// Partition: single line of int32s
func ReadPartitionFromFile(file string) Partition {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	if !sc.Scan() {
		// empty file => empty partition
		return Partition{}
	}
	line := strings.TrimSpace(sc.Text())
	if line == "" {
		return Partition{}
	}

	parts := strings.Fields(line)
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
// Read N and degree list for a CSR (for InitializeCommStats)
func ReadGraphDegreesInput(file string) (int32, []float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadGraphDegreesInput: missing N")
	}
	Nline := strings.TrimSpace(sc.Text())
	N64, err := strconv.ParseInt(Nline, 10, 32)
	if err != nil {
		panic(err)
	}
	N := int32(N64)

	if !sc.Scan() {
		panic("ReadGraphDegreesInput: missing degree line")
	}
	degLine := strings.TrimSpace(sc.Text())
	degParts := strings.Split(degLine, " ")
	degrees := make([]float32, len(degParts))
	for i, p := range degParts {
		fv, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		degrees[i] = float32(fv)
	}
	return N, degrees
}

// Read N, degrees, and partition for InitializeCommStatsFromPartition
func ReadGraphDegreesAndPartitionInput(file string) (int32, []float32, Partition) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadGraphDegreesAndPartitionInput: missing N")
	}
	Nline := strings.TrimSpace(sc.Text())
	N64, err := strconv.ParseInt(Nline, 10, 32)
	if err != nil {
		panic(err)
	}
	N := int32(N64)

	if !sc.Scan() {
		panic("ReadGraphDegreesAndPartitionInput: missing degrees")
	}
	degLine := strings.TrimSpace(sc.Text())
	degParts := strings.Split(degLine, " ")
	degrees := make([]float32, len(degParts))
	for i, p := range degParts {
		fv, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		degrees[i] = float32(fv)
	}

	if !sc.Scan() {
		panic("ReadGraphDegreesAndPartitionInput: missing partition")
	}
	partLine := strings.TrimSpace(sc.Text())
	partParts := strings.Split(partLine, " ")
	P := make(Partition, len(partParts))
	for i, p := range partParts {
		v, err := strconv.ParseInt(p, 10, 32)
		if err != nil {
			panic(err)
		}
		P[i] = int32(v)
	}

	return N, degrees, P
}

// CommStats from 3 lines: Tot, Size, In (or "nil")
func ReadCommStatsFromFile(file string) *CommStats {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadCommStatsFromFile: missing Tot")
	}
	totLine := strings.TrimSpace(sc.Text())
	tot := readFloat32SliceFromLine(totLine)

	if !sc.Scan() {
		panic("ReadCommStatsFromFile: missing Size")
	}
	sizeLine := strings.TrimSpace(sc.Text())
	size := readInt32SliceFromLine(sizeLine)

	if !sc.Scan() {
		panic("ReadCommStatsFromFile: missing In")
	}
	inLine := strings.TrimSpace(sc.Text())
	var in []float32
	if inLine == "nil" || inLine == "" {
		in = make([]float32, len(tot)) // or zeroed slice
	} else {
		in = readFloat32SliceFromLine(inLine)
	}

	return &CommStats{
		Tot:  tot,
		In:   in,
		Size: size,
	}
}

// For EnsureCapacity output (same 3-line format)
func ReadCommStatsPartialFromFile(file string) *CommStats {
	return ReadCommStatsFromFile(file)
}

// For EnsureCapacity input
func ReadEnsureCapacityInput(file string) (*CommStats, int32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadEnsureCapacityInput: missing Tot")
	}
	totLine := strings.TrimSpace(sc.Text())
	tot := readFloat32SliceFromLine(totLine)

	if !sc.Scan() {
		panic("ReadEnsureCapacityInput: missing Size")
	}
	sizeLine := strings.TrimSpace(sc.Text())
	size := readInt32SliceFromLine(sizeLine)

	if !sc.Scan() {
		panic("ReadEnsureCapacityInput: missing In")
	}
	inLine := strings.TrimSpace(sc.Text())
	var in []float32
	if inLine == "nil" || inLine == "" {
		in = nil
	} else {
		in = readFloat32SliceFromLine(inLine)
	}

	if !sc.Scan() {
		panic("ReadEnsureCapacityInput: missing c")
	}
	cLine := strings.TrimSpace(sc.Text())
	c64, err := strconv.ParseInt(cLine, 10, 32)
	if err != nil {
		panic(err)
	}
	c := int32(c64)

	return &CommStats{
		Tot:  tot,
		In:   in,
		Size: size,
	}, c
}

func readFloat32SliceFromLine(line string) []float32 {
	if line == "" {
		return []float32{}
	}
	parts := strings.Split(line, " ")
	out := make([]float32, len(parts))
	for i, p := range parts {
		fv, err := strconv.ParseFloat(p, 32)
		if err != nil {
			panic(err)
		}
		out[i] = float32(fv)
	}
	return out
}

func readInt32SliceFromLine(line string) []int32 {
	if line == "" {
		return []int32{}
	}
	parts := strings.Split(line, " ")
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

// For MoveBuffers.Add and Reset input
func ReadMoveBuffersAddInput(file string) (int, []MoveBuffersOp) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadMoveBuffersAddInput: missing capGuess")
	}
	capLine := strings.TrimSpace(sc.Text())
	capGuess64, err := strconv.ParseInt(capLine, 10, 64)
	if err != nil {
		panic(err)
	}
	capGuess := int(capGuess64)

	if !sc.Scan() {
		panic("ReadMoveBuffersAddInput: missing number of ops")
	}
	kLine := strings.TrimSpace(sc.Text())
	k64, err := strconv.ParseInt(kLine, 10, 64)
	if err != nil {
		panic(err)
	}
	k := int(k64)

	ops := make([]MoveBuffersOp, 0, k)
	for i := 0; i < k; i++ {
		if !sc.Scan() {
			panic("ReadMoveBuffersAddInput: missing op line")
		}
		line := strings.TrimSpace(sc.Text())
		parts := strings.Split(line, " ")
		if len(parts) != 2 {
			panic("ReadMoveBuffersAddInput: op line must have 2 fields")
		}
		c64, err := strconv.ParseInt(parts[0], 10, 32)
		if err != nil {
			panic(err)
		}
		wf, err := strconv.ParseFloat(parts[1], 32)
		if err != nil {
			panic(err)
		}
		ops = append(ops, MoveBuffersOp{
			CommID: int32(c64),
			W:      float32(wf),
		})
	}
	return capGuess, ops
}

// For output of MoveBuffers.Add: two lines in one file.
// line 1: CommIDs
// line 2: CommWts
func ReadInt32SliceFromFile(file string) ([]int32, []float32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadInt32SliceFromFile: missing CommIDs")
	}
	idsLine := strings.TrimSpace(sc.Text())
	ids := readInt32SliceFromLine(idsLine)

	if !sc.Scan() {
		panic("ReadInt32SliceFromFile: missing CommWts")
	}
	wtsLine := strings.TrimSpace(sc.Text())
	wts := readFloat32SliceFromLine(wtsLine)

	return ids, wts
}

// InitializeRefineBuffers input
func ReadRefineBuffersInitInput(file string) (int32, int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadRefineBuffersInitInput: missing N")
	}
	Nline := strings.TrimSpace(sc.Text())
	N64, err := strconv.ParseInt(Nline, 10, 32)
	if err != nil {
		panic(err)
	}
	if !sc.Scan() {
		panic("ReadRefineBuffersInitInput: missing capGuess")
	}
	capLine := strings.TrimSpace(sc.Text())
	cap64, err := strconv.ParseInt(capLine, 10, 64)
	if err != nil {
		panic(err)
	}
	return int32(N64), int(cap64)
}

// RefineBuffers.Reset input
func ReadRefineBuffersResetInput(file string) (int32, int, []int32) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)

	if !sc.Scan() {
		panic("ReadRefineBuffersResetInput: missing N")
	}
	Nline := strings.TrimSpace(sc.Text())
	N64, err := strconv.ParseInt(Nline, 10, 32)
	if err != nil {
		panic(err)
	}
	if !sc.Scan() {
		panic("ReadRefineBuffersResetInput: missing capGuess")
	}
	capLine := strings.TrimSpace(sc.Text())
	cap64, err := strconv.ParseInt(capLine, 10, 64)
	if err != nil {
		panic(err)
	}
	if !sc.Scan() {
		panic("ReadRefineBuffersResetInput: missing queue line")
	}
	qLine := strings.TrimSpace(sc.Text())
	queue := readInt32SliceFromLine(qLine)

	return int32(N64), int(cap64), queue
}

func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}