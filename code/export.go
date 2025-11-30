package main

import (
	"fmt"
	"os"
	"encoding/csv"
	"strconv"
)


func (g *CSR) WriteEdgeListCSV(path string, undirected bool) error {
    f, err := os.Create(path)
    if err != nil {
        return err
    }
    defer f.Close()

    w := csv.NewWriter(f)
    defer w.Flush()

    // header
    if err := w.Write([]string{"src", "dist", "weight"}); err != nil {
        return err
    }

    n := int(g.N)
    for row := 0; row < n; row++ {
        start := int(g.Indptr[row])
        end := int(g.Indptr[row+1])

        for k := start; k < end; k++ {
            col := g.Indices[k]
            wgt := g.Data[k]

            // If graph is undirected + CSR already symmetric,
            // keep only one copy of each edge.
            if undirected {
                // row is an index; col is a NodeID.
                // If NodeID == row index, this dedups;
                // otherwise you can change the condition.
                if int(col) < row {
                    continue
                }
            }

            rec := []string{
                fmt.Sprint(row),              // src index (0-based)
                fmt.Sprint(col),              // dst (NodeID)
                fmt.Sprintf("%g", float64(wgt)), // weight as float
            }
            if err := w.Write(rec); err != nil {
                return err
            }
        }
    }

    return w.Error()
}

func (g *CSR) WriteTripletsCSV(path string, undirected bool) error {
    f, err := os.Create(path)
    if err != nil {
        return err
    }
    defer f.Close()

    w := csv.NewWriter(f)
    defer w.Flush()

    if err := w.Write([]string{"i", "j", "x"}); err != nil {
        return err
    }

    n := int(g.N)
    for row := 0; row < n; row++ {
        start := int(g.Indptr[row])
        end := int(g.Indptr[row+1])

        for k := start; k < end; k++ {
            col := int(g.Indices[k])
            wgt := g.Data[k]

            if undirected && col < row {
                continue
            }

            rec := []string{
                fmt.Sprint(row + 1),      // 1-based
                fmt.Sprint(col + 1),      // 1-based
                fmt.Sprintf("%g", float64(wgt)),
            }
            if err := w.Write(rec); err != nil {
                return err
            }
        }
    }

    return w.Error()
}

func WriteFinalPartitionCSV(sigIds []string, P Partition, path string) error {
    f, err := os.Create(path)
    if err != nil {
        return err
    }
    defer f.Close()

    w := csv.NewWriter(f)
    defer w.Flush()

    // header
    if err := w.Write([]string{"node","node_id", "community"}); err != nil {
        return err
    }

    for i, c := range P {
        rec := []string{
            strconv.Itoa(i),
            sigIds[i],              // node index (must match src/dst in graph_edges.csv)
            strconv.Itoa(int(c)),         // community label
        }
        if err := w.Write(rec); err != nil {
            return err
        }
    }

    return w.Error()
}