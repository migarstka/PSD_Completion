module Completion
    using GraphModule, TreeModule
    export psdCompletion

    # TODO: Check if the order is correct

    function psdCompletion(A)

        # generate graph
        g = Graph(A)
        # find perfect ordering and generate relevant trees to obtain clique tree
        mcsmSearch!(g)
        elimTree = createTreeFromGraph(g)
        superNodeElimTree = createSupernodeEliminationTree(elimTree,g)
        cliqueTree = createCliqueTree(superNodeElimTree,g)
        N = numberOfCliques(cliqueTree)


        # positive semidefinite completion (from Vandenberghe - Chordal Graphs..., p. 362)
        W = copy(A)
        # loop through supernodes in inverse topological  order (order of representative vertex) (i.e. fill W from bottom right to top left)

        for j=(N-1):-1:1

            node = cliqueTree.nodes[cliqueTree.reverseOrder[j]]
            # index set of snd(i) sorted using the numerical ordering i,i+1,...i+ni
            ν = node.value_btm
            # index set containing the elements of col(i) \ snd(i) sorted using numerical ordering σ(i)
            α = node.value_top[end:-1:1]
            # index set containing the the row indizes of the lower-triangular zeros in column i (i: representative index) sorted by σ(i)
            i = node.value_btm[1]
            col_i = node.value_top
            snd = node.value_btm
            η = collect(i+1:1:N)
            #@show(η,α,ν)
            # filter out elements in lower triangular part of column i that are non-zero
            for vertex in col_i
                filter!(f -> f!=vertex,η)
            end
            sort!(η, by=x->g.ordering[x])
            # try factorization first (if matrix has full rank)
            if (rank(W[α,α]) == size(α,1) )
                Y = W[α,α]\W[α,ν]
                W[η,ν] =  W[η,α] * Y
            else
                # use singular value decomposition
                W[η,ν] = W[η,α] * pinv(W[α,α]) * W[α,ν]
            end

            # ensure symmetry
            W[ν,η] = W[η,ν]'
        end

        return W
    end

end # MODULE