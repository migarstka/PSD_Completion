# Test constructors of graph type

workspace()
include("../graph.jl")

using GraphModule, Base.Test

g = Graph([[2, 5, 6], [1, 3, 6], [2, 4, 7], [3, 5, 6, 7], [1, 4], [1, 2, 4], [3, 4]])


mcsmSearch(g)
@test isPerfectOrdering(g) == true



