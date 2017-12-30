# workspace()
# include("../graph.jl")
# include("../tree.jl")

# using GraphModule, TreeModule

# load reference functions from python chompack package
using PyCall, Base.Test
 @pyimport chompack as chom
 @pyimport numpy as np
 @pyimport cvxopt as cvx
 @pyimport scipy.sparse as pysparse
jlmat2pymat(S::SparseMatrixCSC) = pysparse.csc_matrix((S.nzval, S.rowval .- 1, S.colptr .- 1), shape=size(S))
pymat2jlmat(S::PyObject) = SparseMatrixCSC(S[:m], S[:n], S[:indptr] .+ 1, S[:indices] .+ 1, S[:data])



# define example matrix and corresponding graph
A = [
1 -1 2 0 0 1;
-1 2 1 0 0 0;
 2 1 3 0 2 0;
 0 0 0 4 1 0;
 0 0 2 1 5 0;
 1 0 0 0 0 6
];

B = [
  8.0  -0.5   4.0   2.5  0.0   0.0   0.0  0.0   0.0  1.0;
 -0.5   7.5   0.0   5.0  0.0   0.0   5.0  0.0   0.0  0.0;
  4.0   0.0   9.5   4.5  1.5   3.0  -1.0  0.0   0.0  0.5;
  2.5   5.0   4.5  14.0  2.5   0.0   0.0  2.0   2.0  0.0;
  0.0   0.0   1.5   2.5  6.0   0.0   0.0  2.0   0.0  0.0;
  0.0   0.0   3.0   0.0  0.0  10.5  -6.0  0.0   0.0  0.0;
  0.0   5.0  -1.0   0.0  0.0  -6.0  13.0  0.0   0.0  0.0;
  0.0   0.0   0.0   2.0  2.0   0.0   0.0  8.0   0.0  0.0;
  0.0   0.0   0.0   2.0  0.0   0.0   0.0  0.0  10.0  0.0;
  1.0   0.0   0.5   0.0  0.0   0.0   0.0  0.0   0.0  8.0
  ];

  N = size(B,1)
  x = Float64[]
  I = Int64[]
  J = Int64[]

for iii = 1:N
    for jjj = 1:N
        if B[iii,jjj] != 0
            push!(x,B[iii,jjj])
            push!(I,iii-1)
            push!(J,jjj-1)
        end
    end
end

x = [0.522 0.494 0.494 0.5 0.336 0.336 0.4]
#I = [0, 1, 3, 1, 5, 2, 6, 3, 4, 5, 4, 5, 6, 5, 6]
#J = [0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6]
I = [0 1 0 1 2 1 2]
J = [0 0 1 1 1 2 2]
A = cvx.spmatrix(x,(I...),(J...),(3,3))

symb = chom.symbolic(A)
Asymb = chom.cspmatrix(symb)
Asymb += A
W = chom.psdcompletion(Asymb)

(m,n) = W[:size]
F = zeros(m,n)

for i=1:m
    for j=1:n
        F[i,j] = W[i,j]
    end
end

#@test isposdef(F)