# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule, Base.Test

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 1

function validSuperTree(t::Tree,g::Graph)
    return true
end

@testset "Tree Derivation from random Matrices" begin

    for iii=1:nn
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end
        # take random dimension
        dim = rand(rng,2:100,1,1)[1]
        density = rand(rng,0.1:0.1:0.6,1,1)[1]
        # create random sparse matrix
        A = sprand(rng,dim,dim,density)
        A = A+A'
        # create graph from A and make Graph chordal
        g = Graph(A)
        mcsmSearch!(g)
        # create Tree from Graph
        t = createTreeFromGraph(g)
        superNodeElimTree = createSupernodeEliminationTree(t,g)
        @test begin
            res = validSuperTree(superNodeElimTree,g) == true
            if !res
                println(iii)
                println(g)
                println(superNodeElimTree)
            end
            res
        end
    end
end





