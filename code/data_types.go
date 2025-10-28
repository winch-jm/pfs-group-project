package main

type OrderedPair struct{
	x float64
	y float64
}

type Graph struct {
	V []*Node
}
type Matrix [][]float64

type Community struct []*Node

type Partition []*Community 