function testS(S)
  push!(S,1)
  return nothing
end

S= [1,2,2]

testS(S)

println(S)