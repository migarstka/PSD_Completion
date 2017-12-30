workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")

using Helper, Base.Test, GraphModule, TreeModule

rng = MersenneTwister(123554);



# Define number of test matrices
nn = 1000

@testset "Generate Pos-Def-Completable Matrices" begin
    for iii=1:nn

        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.9)
        X = generateCompleatableMatrix(dim,density,rng)

        # check nonzero diagonals
        numDiagZeros = size(find(i->X[i,i] == 0,1:dim),1)
        @test numDiagZeros == 0

        #check transpose
        @test X == X'
    end

end



