workspace()
include("../helper.jl")

using Helper, Base.Test

rng = MersenneTwister(123554);


function generatePosDefMatrix(n::Int64,rng)
    X = rand(rng,n,n)
    X = 1/2*(X+X')
    X = X + n*eye(n)
    return X
end

# Define number of test matrices
nn = 10000

@testset "Test Chordal Extension MCS-M Algorithm" begin

    for iii=1:nn
      dim = rand(rng,2:100,1,1)[1]
      X = generatePosDefMatrix(dim,rng)

      @test isposdef(X) == true
    end

end