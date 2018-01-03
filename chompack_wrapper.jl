# module that wraps the compack psdcompletion to make it easier to use in julia

module chomWrap
    using PyCall
    @pyimport chompack as chom
    @pyimport numpy as np
    @pyimport cvxopt as cvx
    @pyimport scipy.sparse as pysparse
    export completeChompack

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

        symb = chom.symbolic(A, p=cvx.amd[:order])
        Asymb = chom.cspmatrix(symb)
        Asymb += A
        W = chom.psdcompletion(Asymb,reordered=false)
        Z = chom.psdcompletion(Asymb,reordered=true)
        (m,n) = W[:size]
        F = zeros(m,n)
        E = zeros(m,n)

        for i=1:m
            for j=1:n
                F[i,j] = W[i,j]
            end
        end
        for i=1:m
            for j=1:n
                E[i,j] = Z[i,j]
            end
        end
        return F,symb,E
    end


end