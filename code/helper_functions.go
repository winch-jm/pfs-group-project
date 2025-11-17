package main

import (
	// "fmt"
)


func CountCommunities(P Partition) int {
    seen := make(map[int32]struct{})
    for _, c := range P {
        seen[c] = struct{}{}
    }
    return len(seen)
}
