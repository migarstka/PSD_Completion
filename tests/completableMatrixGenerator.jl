workspace()
include("../helper.jl")

using Helper, Base.Test

rng = MersenneTwister(123554);



# Define number of test matrices
nn = 10000

@testset "Generate  Pos-Def-Completable Matrices" begin
    for iii=1:nn

        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.9)
        X,Y = generateCompleatableMatrix(dim,density,rng)

        # check density
        numZeros = size(find(a->a==0,X),1)
        desiredNumZeros = Int(floor((dim^2-dim)*(1-density)))
        desiredNumZeros = desiredNumZeros - mod(desiredNumZeros,2)
        @test numZeros == desiredNumZeros

        # check nonzero diagonals
        numDiagZeros = size(find(i->X[i,i] == 0,1:dim),1)
        @test numDiagZeros == 0

        #check transpose
        @test X == X'
    end

end



