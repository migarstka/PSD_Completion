# Test MCS-M Algorithm with a set of random sparse matrices

workspace()
include("../graph.jl")

using GraphModule, Base.Test

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 3000


@testset "Test Chordal Extension MCS-M Algorithm" begin

    for iii=1:nn
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end
        # take random dimension
        dim = rand(rng,2:100,1,1)[1]
        density = rand(rng,0.1:0.1:0.6,1,1)[1]
        # create random sparse matrix
        A = sprand(rng,dim,dim,density)
        A = A+A'
        # create graph from A
        g = Graph(A)
        # find chordal extension and perfect elimination ordering
        #println(g)
        mcsmSearch!(g)
        @test begin
            res = isPerfectOrdering(g) == true
            if !res
                println(iii)
                println(g)
            end
            res
        end
    end
end





