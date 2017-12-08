# Test constructors of graph type

workspace()
include("../graph.jl")

using GraphModule, Base.Test

g = Graph([[2, 4, 5,7], [1,6,8], [4,5,7], [1,3,6], [1,3,6], [2,4,5], [1,3,8], [2,7]])

mcsmSearch!(g)
@test isPerfectOrdering(g) == true



