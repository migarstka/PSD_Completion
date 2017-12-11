# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")


using GraphModule, TreeModule, Base.Test, Completion, Helper


# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(222);
# Define number of test matrices

# TODO: Change back helper random psd matrix completion range
# take random dimension
nn = 3
@testset "Test positive semidefinite completion" begin

    for iii = 1:nn
        dim = rand(rng,2:100,1,1)[1]
        density = rand(rng,0.1:0.1:0.9,1,1)[1]
        A = generateCompleatableMatrix(30,density,rng)
        @test isposdef(A) == true
    end
end