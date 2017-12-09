# Test constructors of graph type

workspace()
include("../graph.jl")
include("../tree.jl")

using GraphModule, Base.Test, TreeModule


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

#g = Graph([[2, 4, 5,7], [1,6,8], [4,5,7], [1,3,6], [1,3,6], [2,4,5], [1,3,8], [2,7]])

# create random sparse matrix
A = sprand(15,15,0.5)
A = A+A'

# create graph from A and make Graph chordal
g = Graph(A)
mcsmSearch!(g)
# create Tree from Graph
t = createTreeFromGraph(g)
superNodeElimTree = createSupernodeEliminationTree(t,g)
ct = createCliqueTree(superNodeElimTree,g)


checkCliques(ct,g)



