# Test routine to check if the algorithm can complete the randomly generated completable test matrices
# stored in completableMatrices.jld, and created by the createTestMatrices.jl file

workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")

using GraphModule, TreeModule, Base.Test, Completion, Helper, JLD


# load test matrices from file
dict = load("completableMatrices.jld")
testMatrices = dict["randomMatrices"]
numberSaved = size(testMatrices,1)

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

    for iii=1:numberSaved
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end

        # load test matrix
        A = testMatrices[iii]

        # perform positive definite completion with own implementation
        W = psdCompletion(A)

        # check for positive definiteness of result
        @test isposdef(W) == true

        @test matchingMatrices(A,W) == true

    end
end


