# Test constructors of graph type

workspace()
include("../graph.jl")

using GraphModule, Base.Test

g = Graph([[2, 5, 6, 8, 3], [1, 4, 7], [1], [2, 6], [1, 8], [1, 4, 7, 8], [2, 6, 8], [1, 5, 6, 7]])

mcsmSearch!(g)
@test isPerfectOrdering(g) == true



