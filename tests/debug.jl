# Test file to find out why the testMatrix in testMatrix.jld cant be completed by my algorithm
# Ovservations:
# - Only filled 0 and sometimes even overwritten nonzero entries

workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")
include("../completion.jl")
include("../chompack_wrapper.jl")

using GraphModule, Base.Test, TreeModule, Helper, Completion,JLD,chomWrap

dict = load("testMatrix.jld")
A = dict["A"]

W,g = psdCompletion(A)
Q,symb,E = completeChompack(A)
@testset "test" begin
    @test isposdef(W) == true
    @test isposdef(Q) == true
end