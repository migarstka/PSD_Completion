workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")

using GraphModule, TreeModule,Helper, Base.Test

rng = MersenneTwister(222);


# Define number of test matricesw
nn = 1000

@testset "Test Random Positive Definite Symmetric Matrix Generator" begin

    for iii=1:nn
      dim = rand(rng,2:100,1,1)[1]
      X = generatePosDefMatrix(dim,rng)

      @test isposdef(X) == true
      @test X == X'

      # check nonzero diagonals
        numDiagZeros = size(find(i->X[i,i] == 0,1:dim),1)
        @test numDiagZeros == 0
    end

end