# Test routine to check if the algorithm can complete the randomly generated completable test matrices
# genereate on the fly

workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")

using GraphModule, TreeModule, Base.Test, Completion, Helper, JLD


nn = 10000
rng = MersenneTwister(63781);

function matchingMatrices(A,W)
    m,n = size(A)
    result = true
    for iii = 1:m
        for jjj = iii:n
            if A[iii,jjj] != 0.0
                if A[iii,jjj] != W[iii,jjj]
                    result = false
                end
            end
        end
    end
    return result
end

@testset "Test positive semidefinite completion" begin

    for iii=1:nn
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end

        # generate test matrix
        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.1:0.8)
        A = generateCompleatableMatrix(dim,density,rng)


        # perform positive definite completion with own implementation
        W = psdCompletion(A)

        # check for positive definiteness of result
        @test isposdef(W) == true

        @test matchingMatrices(A,W) == true

    end
end


