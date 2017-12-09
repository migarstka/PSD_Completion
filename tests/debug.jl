# Test constructors of graph type

workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")
include("../completion.jl")

using GraphModule, Base.Test, TreeModule, Helper, Completion

rng = MersenneTwister(123554);


#g = Graph([[2, 4, 5,7], [1,6,8], [4,5,7], [1,3,6], [1,3,6], [2,4,5], [1,3,8], [2,7]])

# create random completable matrix
 A, B = generateCompleatableMatrix(8,0.1,rng)
# create graph from A and make Graph chordal
 # generate graph
g = Graph(A)
# find perfect ordering and generate relevant trees to obtain clique tree
mcsmSearch!(g)
elimTree = createTreeFromGraph(g)
superNodeElimTree = createSupernodeEliminationTree(elimTree,g)
cliqueTree = createCliqueTree(superNodeElimTree,g)
N = numberOfCliques(cliqueTree)

W = copy(A)




