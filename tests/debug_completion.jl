# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
#include("../helper.jl")


using PyCall, GraphModule, TreeModule, Base.Test, Completion, Helper


# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(222);
# Define number of test matrices

# TODO: Change back helper random psd matrix completion range
# take random dimension
dim = 20
density = 0.1
A = generateCompleatableMatrix(dim,density,rng)



# create graph from A and make Graph chordal
W = psdCompletion(A)
C = completeChompack(A)
#@test isposdef(C) == true
@test isposdef(W) == true


