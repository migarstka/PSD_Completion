
workspace()
include("../graph.jl")
include("../tree.jl")
include("../completion.jl")
include("../helper.jl")

using GraphModule, TreeModule, Base.Test, Completion, Helper, JLD

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);


# load test matrices from file
dict = load("verifiedMatrizes.jld")
verifiedMatrizes = dict["verifiedMatrizes"]
numberSaved = size(verifiedMatrizes,1)
#testM = [10 4 0 0 0 -4;4 10 -2 -2 0 -2;0 -2 12 8 0 0; 0 -2 8 9 -1 2; 0 0 0 -1 7 -4; -4 -2 0 2 -4 9]
unsuccessful = []
sizes = []

@testset "Test positive semidefinite completion" begin

    for iii=1:numberSaved
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end

        # load test matrix
        A = verifiedMatrizes[iii]

        # create graph from A and make Graph chordal
        W,g = psdCompletion(A)

        @test isposdef(W) == true
        if !isposdef(W)
            push!(unsuccessful,iii)
            push!(sizes,size(A,1))
        end
    end
    # save to JLD file
filename = "testMatrix.jld"
JLD.save(filename, "A", verifiedMatrizes[34])
end
println(unsuccessful)

#save smallest not working test matrix




