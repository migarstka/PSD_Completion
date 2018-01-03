# Test speed of MCS-M Algorithm
workspace()
include("../graph.jl")

using GraphModule, Base.Test

# precompile functions in question
g = Graph([[2, 4, 5,7], [1,6,8], [4,5,7], [1,3,6], [1,3,6], [2,4,5], [1,3,8], [2,7]])
mcsmSearch!(g)


# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(123554);
# Define number of test matrices
nn = 3000

times = zeros(nn)

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
    times[iii] = @elapsed mcsmSearch!(g)
    #println(isPerfectOrdering(g))
end

# determine average runtime for mcsmSearch!
println("Average runtime of mcsmSearch!: $(mean(times)) s")




