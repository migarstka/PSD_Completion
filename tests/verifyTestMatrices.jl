# Test Algorithm to create a superNodeElimTree from a Random Matrix

workspace()
include("../chompack_wrapper.jl")

using Base.Test, chomWrap, JLD

# load test matrices from file
dict = load("completableMatrices.jld")
randomMatrices = dict["randomMatrices"]
numberSaved = dict["numberSaved"]

verifiedMatrizes = []
@testset "Test positive semidefinite completion" begin

    for iii=1:numberSaved
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end

        # load test matrix
        A = randomMatrices[iii]

        # create graph from A and make Graph chordal
        W,symb,E = completeChompack(A)

        #@test isposdef(W) == true
        if isposdef(W)
            push!(verifiedMatrizes,A)
        end
    end
end



# save to JLD file
filename = "verifiedMatrizes.jld"
JLD.save(filename, "verifiedMatrizes", verifiedMatrizes)
