workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule

# create test matrix based on following chordal graph
g = Graph([[2,3,4,5],[1,5],[1,4],[1,3,5,6],[1,2,4],[4]])

# find cliques
mcsSearch!(g)
t = createTreeFromGraph(g)
snet = createSupernodeEliminationTree(t,g)
ct = createCliqueTree(snet,g)

#create matrix
X = zeros(6,6)
for iii = 1:size(ct.nodes,1)
    indizes = ct.nodes[iii].
end