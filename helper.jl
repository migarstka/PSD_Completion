module Helper
    using GraphModule,TreeModule
    export generatePosDefMatrix, generateCompleatableMatrix, completeChompack


    function generatePosDefMatrix(n::Int64,rng)
        X = rand(rng,n,n)
        X = 0.5*(X*X')
        return X
    end

    function droplower(A::SparseMatrixCSC)
        m,n = size(A)
        rows = rowvals(A)
        vals = nonzeros(A)
        V = Vector{eltype(A)}()
        I = Vector{Int}()
        J = Vector{Int}()
        for i=1:n
            for j in nzrange(A,i)
                rows[j]>i && break
                push!(I,rows[j])
                push!(J,i)
                push!(V,vals[j])
            end
        end
        return sparse(I,J,V,m,n)
    end

    """
    generateCompleatableMatrix(n::Int64,density::Float64,rng)

    Generates a n x n-matrix X with specified density that is guaranteed to be completable to a
    positive definite matrix by selecting other values for the zeros. The matrix is created by starting
    with a positive definite matrix Y and then randomly overwriting off-diagonal elements with zeros until
    the desired density is achieved.
    """
    function generateCompleatableMatrix(n::Int64,density::Float64,rng)
        if density <= 0 || density > 1
            error("Density value has to be between 0 and 1")
        end

        # create a random sparse matrix of dimension n and density
        X = sprand(rng,n,n,density) + sparse(diagm(rand(n)))
        X = droplower(X)
        X = full(Symmetric(X))

        # calculate cliques of X
        g = Graph(X)
        mcsmSearch!(g)
        t = createTreeFromGraph(g)
        snet = createSupernodeEliminationTree(t,g)
        ct = createCliqueTree(snet,g)
        println(ct)
        # overwrite the cliques with positive definite matrices and set the rest to zero
        Y = 0*X
        for clique in ct.nodes
            vertices = clique.value_btm
            N = size(vertices,1)
            Y[vertices,vertices] = generatePosDefMatrix(N,rng)
        end

        return Y

    end

    function completeChompack(A)

        N = size(A,1)
        x = Float64[]
        I = Int64[]
        J = Int64[]

        for iii = 1:N
            for jjj = 1:N
                if A[iii,jjj] != 0
                    push!(x,A[iii,jjj])
                    push!(I,iii-1)
                    push!(J,jjj-1)
                end
            end
        end

        A = cvx.spmatrix(x,(I...),(J...))

        symb = chom.symbolic(A)
        Asymb = chom.cspmatrix(symb)
        W = chom.psdcompletion(Asymb)

        (m,n) = W[:size]
        F = zeros(m,n)

        for i=1:m
            for j=1:n
                F[i,j] = W[i,j]
            end
        end
        return F
    end
end #Module
