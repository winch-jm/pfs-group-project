package main

// unique int id for cell (node)
type NodeID = int32      // up to ~2B cells

// used to index CSR arrays, marks where adjacency list starts and ends
type Idx    = int32
// edge weight
type Weight = float32    // float32 is plenty for similarities/distances


// Preprocessing only (computing similarities / build KNN --> discard)
type DenseRows struct {
    // n rows (cells) × d columns (genes or features)
    N, D int
    Data []float32 // row-major: row i at Data[i*D:(i+1)*D]
}

// Compressed Sparse Row graph (undirected; store symmetric edges)
// Used from the first KNN graph onward, at every level
// - fast neighbor scans in local-move & refinement
type CSR struct {
    N      int32     // #nodes
    Indptr []Idx     // len N+1; row i neighbors => Indices[Indptr[i]:Indptr[i+1]]
	Indices []NodeID  // neighbor ids
    Data   []Weight  // edge weights parallel to Indices
    // Cached degrees (sum of weights) – required by modularity/CPM bookkeeping
    Degree []float32 // len N;  supports deltaQ math
    TwoM float32
}

// Used throughout, P[i] = community id of node i;
// mutated during local move & refinement; remapped after aggregation
type Partition []int32


// Used during local move for deltaQ bookkeeping 
// updated when node leaves/enters community
type CommStats struct {
    // Sum of degrees of nodes in each community (2m_c in undirected weighted graphs)
    Tot []float32   // len = current #communities

    // (Optional) Internal weight per community if you use quality functions needing it
    In  []float32   // len = current #communities

	// #nodes per community
	Size []int32   // len = current #communities
}

// Used for deltaQ calculations
type QualityFn int
const (
    Modularity QualityFn = iota
    CPM                     // Constant Potts Model
	RBPM					// Reichardt Bornhold Potts Model
)
// Used for deltaQ calculations
type QualityCfg struct {
    Fn        QualityFn
    Gamma     float32  // resolution parameter (γ); 1.0 for classic modularity
    SelfLoop  float32  // small self-loop to stabilize (e.g., 1e-3 or avg weight / k)
}

// Used per node during local-move 
// when calculating 
type MoveBuffers struct {
    // Map neighbor community -> total edge weight from current node
    // Reused per node to avoid allocations.
    CommIDs   []int32 // unique community ids seen around node i
    CommWts   []float32 // parallel; sum of weights to each community (k_i,in(c))
    Seen      map[int32]int // community id -> index in CommIDs/CommWts

	Order   []int // indices into CommIDs for argmax
}

// Used when each level completes, stats for stopping criteria and quality metrics
type LevelStats struct {
    Level          int
    Quality        float64
    NumCommunities int
    Moves          int
}

// map old node to supernode id after refinement, to contract graph and remap partition
type AggMap []int32 // len = N; AggMap[i] = supernode id of node i after refinement

// Coordinate List
// Used to build next-level graph from inter-community edges
type COO struct {
    Row []int32
    Col []int32
    Val []float32
}

// Used for refinement step (BFS/DFS to split weakly/poorly connected parts)
type RefineBuffers struct {
    Queue   []int32   // BFS/DFS queue
    Visited []int32   // visitation stamp per node (avoid allocating bool)
    Stamp   int32
}

type IndexScore struct{
    similarity float32
    index int32
    heapIdx int
}