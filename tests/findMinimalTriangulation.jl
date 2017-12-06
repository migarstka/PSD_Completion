
# Test the MCS-M algorithm to find minimal triangulation of the non-chordal input graph
workspace()
include("../graph.jl")

using GraphModule
using Base.Test

# graph from MCS-M paper
println("Graph 1: MCS-M paper\n")
g1 = Graph([[2,10],[1,3,4],[2,5,9],[2,5],[4,6,9],[5,7],[6,8],[5,7,9],[3,5,8,10],[1,9]])
# try to find perfect ordering
mcsSearch(g1)
# check if ordering is perfect, i.e. graph is chordal
println("Graph is chordal: $(isPerfectOrdering(g1))\n")
# fill graph
mcsmSearch(g1)

# circle graph
println("Graph 2: Circular graph\n")
g2 = Graph([[2,8],[1,3],[2,4],[3,5],[4,6],[5,7],[6,8],[1,7]])
mcsSearch(g2)
println("Graph is chordal: $(isPerfectOrdering(g2))\n")
mcsmSearch(g2)


# star graph
println("Graph 3: Star graph\n")
g3 = Graph([[6],[6],[6],[6],[6],[1,2,3,4,5]])
mcsSearch(g3)
println("Graph is chordal: $(isPerfectOrdering(g3))\n")
mcsmSearch(g3)

# christmas houses
println("Graph 4: Christmas Houses\n")
g4 = Graph([[2,11],[1,3],[2,4],[3,5,11],[4,6],[5,7,10],[6,8],[7,9],[8,10],[6,9,11],[1,4,10]])
mcsSearch(g4)
println("Graph is chordal: $(isPerfectOrdering(g4))\n")
mcsmSearch(g4)

# # Perform Test
@testset "Find Minimal Triangulation Test" begin
    @test isPerfectOrdering(g1) == true
    @test isPerfectOrdering(g2) == true
    @test isPerfectOrdering(g3) == true
    @test isPerfectOrdering(g4) == true
end