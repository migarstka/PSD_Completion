# Test base types and functions

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule


# create graph from adjacency list and test functions
#g = Graph([[2,3,5],[1,3,4,5],[1,2,4,5,7],[2,3,5,6],[1,2,3,4,6,7,8],[4,5],[3,5,8,9],[5,7,9],[7,8]])
g = Graph([[2,3,4,5],[1,5],[1,4],[1,3,5,6],[1,2,4],[4]])
mcsSearch!(g)
t = createTreeFromGraph(g)
snet = createSupernodeEliminationTree(t,g)
ct = createCliqueTree(snet,g)
