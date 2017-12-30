module Completion
    using GraphModule, TreeModule
    export psdCompletion

    # TODO: Check if the order is correct
    # - Apply Permutation A[p,p] reordering of rows and columns of A
    # - complete matrix and reverse the reordering operation again
    function psdCompletion(A)
        doPrint = true
        # generate graph
        g = Graph(A)
        # find perfect ordering and generate relevant trees to obtain clique tree
        mcsmSearch!(g)
        elimTree = createTreeFromGraph(g)
        superNodeElimTree = createSupernodeEliminationTree(elimTree,g)
        cliqueTree = createCliqueTree(superNodeElimTree,g)
        Nc = numberOfCliques(cliqueTree)
        N = numberOfVertices(g)
        # get permutation vector p and compute inverse permutation vector ip
        p = g.ordering
        ip = zeros(Int64,length(p))
        ip[p] = 1:length(p)

        # positive semidefinite completion (from Vandenberghe - Chordal Graphs..., p. 362)
        # permutate matrix based on ordering p (p must be a vector type)

        W = A[p,p]
        # loop through supernodes in inverse topological  order (order of representative vertex) (i.e. fill W from bottom right to top left)

        for j=(Nc-1):-1:1
            doPrint && println(">>> j = $(j) current W:")
            doPrint && @show(W)
            node = cliqueTree.nodes[cliqueTree.reverseOrder[j]]
            # index set of snd(i) sorted using the numerical ordering i,i+1,...i+ni
            ν = g.ordering[node.value_btm]
            # index set containing the elements of col(i) \ snd(i) sorted using numerical ordering σ(i)
            α = g.ordering[node.value_top]
            # index set containing the row indizes of the lower-triangular zeros in column i (i: representative index) sorted by σ(i)
            i = ν[end]
            η = collect(i+1:1:N)
            doPrint && println("Node=$(node)\nν=$(ν)\nα=$(α)\ni=$(i)\nη_before=$(η)")

            # filter out elements in lower triangular part of column i that are non-zero
            filter!(x -> !in(x,α),η)

            doPrint && println("η_filtered = $(η)")
            # FIXME: Should be already appropriately sorted
            #sort!(η, by=x->g.ordering[x])
            # try factorization first (if matrix has full rank)
            Waa = W[α,α]
            if (rank(Waa) == size(α,1) )

                doPrint && println("Full Rank factorization with W[α,α]=$(W[α,α]), W[α,ν]=$(W[α,ν])")
                Y = Waa\W[α,ν]
                W[η,ν] =  W[η,α] * Y
            else
                doPrint && println("SVD factorization with W[α,α]=$(W[α,α]), W[α,ν]=$(W[α,ν])")

                # otherwise use eigenvalue decomposition to compute pseudo inverse
                F = eigfact(Waa)
                Z = F[:vectors]
                λ = F[:values]
                tol = 1e-15 * max(λ)
                for iii in size(λ,1)
                    ev = λ[iii]
                    if ev < tol
                        λ[iii] = 0.0
                    else
                        λ[iii] = 1.0/ev
                    end
                end
                Waa_pinv = Z*diagm(λ)*Z'
                W[η,ν] = W[η,α] * Waa_pinv * W[α,ν]
            end

            # symmetry condition
            W[ν,η] = W[η,ν]'
        end

        # invert the permutation
        W = W[ip,ip]

        return W,g,superNodeElimTree, cliqueTree
    end

end # MODULE