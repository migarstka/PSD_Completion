# Test Algorithm to check the completable matrices with the chompack library
# As it turns out only some of the matrices from completableMatrices.jld can be
# completed by the chompack routine
# -> Those are stored in a verifiedMatrizes array, the other cases are stored in  an unverifiedMatrizes array
# and exported using JLD

workspace()
include("../chompack_wrapper.jl")

using Base.Test, chomWrap, JLD

# load test matrices from file
dict = load("completableMatrices.jld")
randomMatrices = dict["randomMatrices"]
numberSaved = dict["numberSaved"]

verifiedMatrizes = []
unverifiedMatrizes = []
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
        else
            push!(unverifiedMatrizes,A)
        end
    end
end



# save to JLD file
filename = "verifiedMatrizes.jld"
JLD.save(filename, "verifiedMatrizes", verifiedMatrizes,"unverifiedMatrizes",unverifiedMatrizes)
