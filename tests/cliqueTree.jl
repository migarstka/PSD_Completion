# Test Algorithm to create a clique tree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, TreeModule, Base.Test

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 1000


function checkLowerRow(ct::Tree,superTree::Tree)
    # the value_btm property of the clique tree has to coincide with the value_top property of the superNode elimination tree
    for iii=1:size(ct.nodes,1)
        if ct.nodes[iii].value_btm != superTree.nodes[iii].value_top
            return false
        end
    end
    return true
end

function checkUpperRow(ct::Tree,superTree::Tree,g::Graph)
    for node in ct.nodes
        if node.parent != 0
            snd = node.value_btm
            col = Int64[]
            for v in snd
                # find higher degree neighbors of each v
                neighbors = g.adjacencyList[v]
                for neighbor in neighbors
                    if g.ordering[neighbor] > g.ordering[v]
                        push!(col,neighbor)
                    end
                end
            end
            # only consider the unique values
            col = union(col)
            # remove values from snd
            col_snd = filter(x->!in(x,snd),col)

            if size(node.value_top,1) != size(col_snd,1)
                return false
            end

            for v in node.value_top
                if !in(v,col_snd)
                    return false
                end
            end
        end
    end
    return true
end

# all vertices in a clique have to be pair-wise adjacent, i.e. induce a complete subgraph
function checkCliques(ct::Tree,g::Graph)
    for node in ct.nodes
        clique = node.value_btm
        if size(clique,1) > 1
            for v in clique
                otherNodes = filter(x->x!=v,clique)
                for cliqueNode in otherNodes
                    if !in(cliqueNode,g.adjacencyList[v])
                        return false
                    end
                end
            end
        end
    end
    return true
end


@testset "Clique Tree Derivation from random Matrices" begin

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
        ct = createCliqueTree(snet,g)
        @test checkLowerRow(ct,snet)
        @test   checkUpperRow(ct,snet,g)
        @test checkCliques(ct,g)

    end
end





