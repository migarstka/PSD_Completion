 workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")
include("../chompack_wrapper.jl")
 include("../completion.jl")

using Helper, chomWrap, Base.Test,Completion

rng = MersenneTwister(123554);

X,g = generateCompleatableMatrix(10,0.5,rng)

W,symb, Z = completeChompack(X)
Q,g_ret,snet, ct = psdCompletion(X)

@testset "Test both" begin
    @test isposdef(W) == true
    @test isposdef(Q) == true
end
# symb[:cliques]()
# symb[:supernodes]()