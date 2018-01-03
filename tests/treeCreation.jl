# Test MCS-M Algorithm with a set of random sparse matrices

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule, Base.Test

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 1000

function validTree(t::Tree,g::Graph)
    N = numberOfVertices(g)
    if size(t.nodes,1) != size(g.adjacencyList,1)
        warn("Number of tree nodes and graph Vertices doesnt match!")
        return false
    end

    if t.root != g.reverseOrder[N]
        warn("Root is not chosen correctly!")
        return false
    end

    # check order of parents of all non-root nodes
    for node in t.nodes
        if node.value_top != t.root
            if g.ordering[node.parent] <= g.ordering[node.value_top]
                warn("Test failed at node $(node)")
                return false
            end
        end

    end
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

        @test begin
            res = validTree(t,g) == true
            if !res
                println(iii)
                println(g)
                println(t)
            end
            res
        end
    end
end





