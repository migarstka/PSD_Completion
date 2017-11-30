
# Test the MCS-M algorithm to find minimal triangulation of the non-chordal input graph
workspace()
include("../graph.jl")

using GraphModule
using Base.Test

# graph from MCS-M paper
g = Graph([[2,10],[1,3,4],[2,5,9],[2,5],[4,6,9],[5,7],[6,8],[5,7,9],[3,5,8,10],[1,9]])
mcsmSearch(g)


# # example with chordal graph -> mcs should lead to perfect elimination ordering
# g_chord = Graph([[2,3,5],[1,3,4,5],[1,2,4,5,7],[2,3,5,6],[1,2,3,4,6,7,8],[4,5],[3,5,8,9],[5,7,9],[7,8]])
# mcsSearch(g_chord)
# g_chord_check = isPerfectOrdering(g_chord)


# # example of a non-chordal graph (edge between 7 and 8 removed)
# g = Graph([[2,3,5],[1,3,4,5],[1,2,4,5,7],[2,3,5,6],[1,2,3,4,6,7,8],[4,5],[3,5,9],[5,9],[7,8]])
# mcsSearch(g)
# g_check = isPerfectOrdering(g)

# # example of a non-chordal graph with very few edges (edge between 7 and 8 removed)
# g_few = Graph([[2],[1,4],[5,7],[2,5],[3,4,6,8],[5],[3,9],[5,9],[7,8]])
# mcsSearch(g_few)
# g_few_check = isPerfectOrdering(g_few)


# # Perform Test
# @testset "Find Minimal Triangulation Test" begin
#     @test g_chord_check == true
#     @test g_check == false
#     @test g_few_check == false

# end