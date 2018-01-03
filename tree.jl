module TreeModule
    using GraphModule
    export Tree, Node, createTreeFromGraph, createSupernodeEliminationTree, createCliqueTree, numberOfCliques


    # -------------------------------------
    # TYPE DEFINITIONS
    # -------------------------------------
    type Node
        value_top::Array{Int64,1}
        value_btm::Array{Int64,1}
        degree::Int64
        parent::Int64
        children::Array{Int64}
        inSuperNode::Int64

        # two constructor definitions
        function Node(value_top::Array{Int64},value_btm::Array{Int64},degree::Int64,parent::Int64)
            new(value_top,value_btm,degree,parent,[],0)
        end

        function Node(value_top::Array{Int64},value_btm::Array{Int64},degree::Int64,parent::Int64,children::Array{Int64})
            new(value_top,value_btm,degree,parent,children,0)
        end

    end

    type Tree
        root::Int64
        nodes::Array{Node}
        order::Array{Int64}
        reverseOrder::Array{Int64}
        #constructor
        function Tree()
            new(0,Int64[],Int64[],Int64[])
        end
    end

# Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::Tree)
    println(io,"\nTree - Nodes:\nRoot: $(obj.root)\nOrder: $(obj.order)\n reverseOrder: $(obj.reverseOrder)")
    for node in obj.nodes
        println("Node - Value_top: $(node.value_top), Value_btm: $(node.value_btm), Degree: $(node.degree), Parent: $(node.parent), Children: $(node.children), inSuperNode: $(node.inSuperNode)\n")
    end
  end


    # -------------------------------------
    # FUNCTION DEFINITIONS
    # -------------------------------------

     function numberOfCliques(ct::Tree)
        return size(ct.nodes,1)
    end


    function createTreeFromGraph(g::Graph)
        tree = Tree()
        N = numberOfVertices(g)
        # loop over Vertices of graph
        for i=1:N
            value = i
            # number of i-neighbors with order higher than order of node i
            higherNeighbors = findHigherNeighbors(g,i)
            degree = size(higherNeighbors,1)
            parent =findParent(g,higherNeighbors)
            order = g.ordering[i]
            node = Node([value],Int64[],degree,parent)
            push!(tree.nodes,node)
            if parent == 0
                tree.root = i
            end
        end

        # fill the children property of each node
        for node in tree.nodes
            if node.parent != 0
                push!(tree.nodes[node.parent].children,node.value_top[1])
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
                # if vertex is representative, i.e. doesnt have lower degree child, create new SuperNode
                superNode = Node([nodeInd],[0],-1,-1)
                push!(superTree.nodes,superNode)
                node.inSuperNode = size(superTree.nodes,1)
            else
               # if node is not representative, add to existing supernode that contains that child
                superNode = superTree.nodes[child.inSuperNode]
                node.inSuperNode = child.inSuperNode
                push!(superNode.value_top,nodeInd)

            end
        end
        # determine parent / children relationship between supernodes
        for iii=1:size(superTree.nodes,1)
            superNode = superTree.nodes[iii]
            highestNodeInd = superNode.value_top[indmax(g.ordering[superNode.value_top])]
            highestNode = t.nodes[highestNodeInd]
            if (highestNode.parent == 0)
                superNode.parent = 0
                superTree.root = iii
            else
                superNode.parent = t.nodes[highestNode.parent].inSuperNode
            end
        end
        # fill the children property of each node
        for iii=1:size(superTree.nodes,1)
            superNode = superTree.nodes[iii]
            if superNode.parent != 0
                push!(superTree.nodes[superNode.parent].children,iii)
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
            # snd_v: the vertices of snd_v are ordered consecutively
            snd_v = copy(superNode.value_top)
            col_v_snd_v = Int64[]
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
            sort!(col_v_snd_v, by=x->g.ordering[x])
            node = Node(col_v_snd_v,snd_v,-1,superNode.parent,superNode.children)
            push!(cliqueTree.nodes,node)
        end

        # Order the supernodes in descending order for algorithm to work
        cliqueTree.reverseOrder = collect(1:numberOfCliques(cliqueTree))
        cliqueTree.order = collect(numberOfCliques(cliqueTree):-1:1)

        return cliqueTree
    end


end #MODULE