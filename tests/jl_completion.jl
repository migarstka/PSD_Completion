# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")

using GraphModule, TreeModule, Base.Test, Completion, Helper

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 1

@testset "Test positive semidefinite completion" begin

    for iii=1:nn
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end

        # take random dimension
        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.1:0.6)
        A = generateCompleatableMatrix(8,density,rng)
          if mod(iii,10) == 0
            #println(A)
        end

        # create graph from A and make Graph chordal
        W = psdCompletion(A)

        @test isposdef(W) == true

    end
end


