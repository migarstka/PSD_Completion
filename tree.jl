module TreeModule
    using GraphModule
    export Tree, Node, createTreeFromGraph, createSupernodeEliminationTree, createCliqueTree


    # -------------------------------------
    # TYPE DEFINITIONS
    # -------------------------------------
    type Node
        value_top
        value_btm
        degree::Int64
        parent::Int64
        children::Array{Int64}
        inSuperNode::Int64

        # two constructor definitions
        function Node(value_top,value_btm,degree::Int64,parent::Int64)
            new(value_top,value_btm,degree,parent,[],0)
        end

        function Node(value_top,value_btm,degree::Int64,parent::Int64,children::Array{Int64})
            new(value_top,value_btm,degree,parent,children,0)
        end

    end

    type Tree
        root::Int64
        nodes::Array{Node}
         #constructor
        function Tree()
            new(0,[])
        end
    end

# Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::Tree)
    println(io,"\nTree - Nodes:\nRoot: $(obj.root)\n")
    for node in obj.nodes
        println("Node - Value_top: $(node.value_top), Value_btm: $(node.value_btm), Degree: $(node.degree), Parent: $(node.parent), Children: $(node.children), inSuperNode: $(node.inSuperNode)\n")
    end
  end


    # -------------------------------------
    # FUNCTION DEFINITIONS
    # -------------------------------------
    function createTreeFromGraph(g::Graph)
        tree = Tree()
        N = numberOfVertizes(g)
        # loop over vertizes of graph
        for i=1:N
            value = i
            # number of i-neighbors with order higher than order of node i
            higherNeighbors = findHigherNeighbors(g,i)
            degree = size(higherNeighbors,1)
            parent =findParent(g,higherNeighbors)
            order = g.ordering[i]
            node = Node(value,0,degree,parent)
            push!(tree.nodes,node)
            if parent == 0
                tree.root = i
            end
        end

        # fill the children property of each node
        for node in tree.nodes
            if node.parent != 0
                push!(tree.nodes[node.parent].children,node.value_top)
            end
        end
        return tree

    end

    function createSupernodeEliminationTree(t::Tree,g::Graph)
        superTree = Tree()
        # go through nodes of tree in topological order
        for nodeInd in g.reverseOrder
            node = t.nodes[nodeInd]
            # check if node is representative node (lowest vertex in clique)
            child = hasLowerDegChild(node,t)
            if child == -1
                # if so create new SuperNode
                superNode = Node([nodeInd],NaN,-1,-1)
                push!(superTree.nodes,superNode)
                node.inSuperNode = size(superTree.nodes,1)

            else
                # if not add to existing supernode that contains that child
                superNode = superTree.nodes[child.inSuperNode]
                node.inSuperNode = child.inSuperNode
                push!(superNode.value_top,nodeInd)
            end
        end
        # determine parent / children relationship between supernodes
        for superNode in superTree.nodes
            highestNodeInd = superNode.value_top[indmax(g.ordering[superNode.value_top])]
            highestNode = t.nodes[highestNodeInd]
            if (highestNode.parent == 0)
                superNode.parent = 0
            else
                superNode.parent = t.nodes[highestNode.parent].inSuperNode
            end
        end
        # fill the children property of each node
        i = 0
        for superNode in superTree.nodes
            i+=1
            if superNode.parent != 0
                push!(superTree.nodes[superNode.parent].children,i)
            end
        end

        return superTree
    end

    # checks if node/vertex is representative (lowest degree in clique)
    function hasLowerDegChild(n::Node,t::Tree)
        for childInd in n.children
            child = t.nodes[childInd]
            if !(child.degree < n.degree + 1)
                return child
            end
        end
        return -1
    end
    # FIXME: The roots
    # takes as input a SupernodeEliminationTree and turns it into a clique tree (s. Vandenberghe, Chordal Graphs and SDO, p.287)
    function createCliqueTree(t::Tree,g::Graph)
        cliqueTree = Tree()
        for superNode in t.nodes
            snd_v = superNode.value_top
            col_v_snd_v = []
            for node in snd_v
                higherNeighbors = findHigherNeighbors(g,node)
                if size(higherNeighbors,1) > 0
                    col_v_snd_v = union(col_v_snd_v,higherNeighbors)
                end
            end
            # filter out duplicates and elements of snd(v)
            for vertex in snd_v
                filter!(f -> f!=vertex,col_v_snd_v)
            end

            node = Node(col_v_snd_v,snd_v,-1,superNode.parent,superNode.children)
            push!(cliqueTree.nodes,node)
        end
        return cliqueTree
    end

end #MODULE