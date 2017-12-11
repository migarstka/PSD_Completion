# Test Algorithm to random psd-completable non-posdef test matrices

workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")


using GraphModule,TreeModule, Helper
import JLD

# create random seed (to get reproducable sequence of random numbers)
rng = MersenneTwister(222);

# Define number of test matrices
nn = 1000

# create variable to store non-pos-def matrices
randomMatrices = Array{Float64,2}[]
numberSaved = 0
for iii = 1:nn
    dim = rand(rng,2:100,1,1)[1]
    density = rand(rng,0.1:0.1:0.4,1,1)[1]
    A = generateCompleatableMatrix(dim,density,rng)
    # if the randomly generated compleatable matrix is not posetive definite, save it
    if !isposdef(A)
      push!(randomMatrices,A)
      numberSaved+=1
    end
end

# save to JLD file
filename = "completableMatrices.jld"
JLD.save(filename, "randomMatrices", randomMatrices, "numberSaved", numberSaved)
println("\nRandom Matrix Generation Completed. $(numberSaved) random sparse non-posdef and compleatable matrices are saved in $(filename)\n")
