# workspace()
# include("../graph.jl")
# include("../tree.jl")

# using GraphModule, TreeModule

# load reference functions from python chompack package
# using PyCall
#  @pyimport osqp as osqp
#  @pyimport chompack as chom
#  @pyimport numpy as np
#  @pyimport cvxopt as cvx

#  @pyimport scipy.sparse as pysparse
# jlmat2pymat(S::SparseMatrixCSC) = pysparse.csc_matrix((S.nzval, S.rowval .- 1, S.colptr .- 1), shape=size(S))
# pymat2jlmat(S::PyObject) = SparseMatrixCSC(S[:m], S[:n], S[:indptr] .+ 1, S[:indices] .+ 1, S[:data])



# create random positive definite matrizes
sizeArr = [2^n for n in 1:8] #2,...,256
sparsity = 0.25

counter = 1
testMatrizes = [];
while counter <= 8
    N = sizeArr[counter]
    M = rand(-1:0.1:1,N,N)
    if rank(M) == N

        # create pos-def matrix
        M = M*M'
        # remove blocks from matrix to make it sparse
        numberOfZeroBlocks = round(N^2*(1-sparsity))




        push!(testMatrizes,M)
        counter+=1


    end

    # make
end

# define example matrix and corresponding graph
A = [
1 -1 2 0 0 1;
-1 2 1 0 0 0;
 2 1 3 0 2 0;
 0 0 0 4 1 0;
 0 0 2 1 5 0;
 1 0 0 0 0 6
];

I = [0, 1, 3, 1, 5, 2, 6, 3, 4, 5, 4, 5, 6, 5, 6]
J = [0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6]
A = cvx.spmatrix(0.1,(I...),(J...))

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
