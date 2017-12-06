# Test of  the isPerfectOrdering() function
workspace()
include("../graph.jl")

using GraphModule
using Base.Test

# example with chordal graph -> mcs should lead to perfect elimination ordering
g_chord = Graph([[2,3,5],[1,3,4,5],[1,2,4,5,7],[2,3,5,6],[1,2,3,4,6,7,8],[4,5],[3,5,8,9],[5,7,9],[7,8]])
mcsSearch(g_chord)
g_chord_check = isPerfectOrdering(g_chord)


# example of a non-chordal graph (edge between 7 and 8 removed)
g = Graph([[2,3,5],[1,3,4,5],[1,2,4,5,7],[2,3,5,6],[1,2,3,4,6,7,8],[4,5],[3,5,9],[5,9],[7,8]])
mcsSearch(g)
g_check = isPerfectOrdering(g)

# example of a non-chordal graph with very few edges (edge between 7 and 8 removed)
g_few = Graph([[2],[1,4],[5,7],[2,5],[3,4,6,8],[5],[3,9],[5,9],[7,8]])
mcsSearch(g_few)
g_few_check = isPerfectOrdering(g_few)

# star graph
g_star = Graph([[2,3,4],[1,5,3],[1,4,5,2],[1,3,5],[2,3,4]])


# graph
g2 = Graph([[3],[3,4],[1,4],[2,3,5],[4]])
# Perform Test
@testset "isPerfect Ordering Test" begin
    @test g_chord_check == true
    @test g_check == false
    @test g_few_check == false
    @test g_star == true

end