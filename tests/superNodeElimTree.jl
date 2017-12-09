# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule, Base.Test

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 10000


#= test validity of super Tree
- all Vertices are covered
- no node has more than one parent
- supdernode root node has to contain the root of the elimination tree
- the parent of each supernode is the supernode that contains the parent of the highest order vertex
=#
function checkVertices(t::Tree,g::Graph)
    # all vertizes of original graph have to occur in the supernodes of the super node elimination tree
    N = numberOfVertices(g)
    vertices = collect(1:N)
    for  superNode in t.nodes
        for  node in superNode.value_top
            deleteat!(vertices,findfirst(vertices,node))
        end
    end
    return size(vertices,1) == 0
end

function checkRoot(t::Tree,superTree::Tree)
    # the root of the elimination tree has to be contained in the root of the super node elimination tree
    return in(t.root,superTree.nodes[superTree.root].value_top)
end

function checkParents(t::Tree,superTree::Tree,g::Graph)
    # no supernode has more than one parent
    for superNode in superTree.nodes
        if size(superNode.parent,1) > 1
            return false
        end
    end
    # the parent of each supernode is the supernode that contains the parent of the highest order vertex
    for superNode in superTree.nodes
        if superNode.parent != 0
            highestVertex = superNode.value_top[indmax(g.ordering[superNode.value_top])]
            parentOfVertex = t.nodes[highestVertex].parent
            if !in(parentOfVertex,superTree.nodes[superNode.parent].value_top)
                return false
            end
        end
    end
    return true
end


@testset "SuperNode Tree Derivation from random Matrices" begin

    for iii=1:nn

        # take random dimension
        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.1:0.6)
        # create random sparse matrix
        A = sprand(rng,dim,dim,density)
        A = A+A'

        # create graph from A and make Graph chordal
        g = Graph(A)
        mcsmSearch!(g)
        # create Tree from Graph
        t = createTreeFromGraph(g)
        snet = createSupernodeEliminationTree(t,g)
        @test checkVertices(snet,g)
        @test   checkRoot(t,snet)
        @test checkParents(t,snet,g)

    end
end





