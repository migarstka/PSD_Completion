

using JLD, Base.Test

# create random seed (to get reproducable sequence of random numbers)


# load test matrices from file
dict = load("completableMatrices.jld")
randomMatrices = dict["randomMatrices"]
numberSaved = dict["numberSaved"]

for i = 1:numberSaved
  A = randomMatrices[i]
  println(A == A')
end