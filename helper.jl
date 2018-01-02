module Helper
    using GraphModule,TreeModule
    export generatePosDefMatrix, generateCompleatableMatrix


    function generatePosDefMatrix(n::Int64,rng)
        # FIXME: Remove the range
        #X = rand(rng,-2.0:1.0:2.0,n,n)
        X = rand(rng,n,n)
        # make sure the diagonals are unequal to zero
        for iii = 1:n
            if X[iii,iii] == 0.0
                X[iii,iii] = rand(rng,0.1:1e-2:1)
            end
        end
        #X = 0.5*(X*X')
        X = (X*X')

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
    positive definite matrix by selecting other values for the zeros. Note: The density is only achieved
    roughly since some zeros are filled by overwritting cliques with positive definite blocks.
    """
    function generateCompleatableMatrix(n::Int64,density::Float64,rng)
        if density <= 0 || density > 1
            error("Density value has to be between 0 and 1")
        end


        # TODO:
        # - generate a full positive definite matrix
        # - create a chordal graph with matching density
        # - delete entries that are not considered in the graph
        # - resulting matrix should be positive definite completable

        # create a full symmetric positive definite matrix of size n
        X = generatePosDefMatrix(n,rng)

        # generate a random chordal graph that matches the demanded density
        A = sprand(rng,n,n,density)
        A = droplower(A)
        A = full(Symmetric(A))
        g = Graph(A)
        mcsmSearch!(g)
        # FIXME: Remove the predefined graph structure
        #g = Graph([[2,3],[1,3],[1,2,4,5],[3,5],[3,4,6],[5]])
        # set entries of edges not in G(V,E) to zero
        for i = 2:n
            for j = 1:i-1
                if i != j
                    if !in(j,g.adjacencyList[i])
                        X[i,j] = 0
                        X[j,i] = 0
                    end
                end
            end
        end

        return X,g

    end

end #Module
