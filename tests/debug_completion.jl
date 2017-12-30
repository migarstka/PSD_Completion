# Test Algorithm to create a superNodeElimTree from a Random Matrix

#workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")


using GraphModule, TreeModule, Base.Test, Completion, Helper


# create random seed (to get reproducable sequence of random numbers)
X = [0.522 0.494 0; 0.494 0.5 0.336; 0 0.336 0.4]
W = psdCompletion(X)
@test isposdef(W) == true
