# Test constructors of graph type

workspace()
include("../graph.jl")

using GraphModule, Base.Test


# function that iterates through each non-zero entry of input matrix A and checks if the corresponding indizes appear as edges in the graph g
function checkGraphNodes(g::Graph,A)
    A = full(A)
    for i = 1:size(A,1)
        for j=1:size(A,2)
            if A[i,j] != 0 && i != j
                if !in(j,g.adjacencyList[i]) || !in(i,g.adjacencyList[j])
                    return false
                end
            end
        end
    end
    return true
end

# Define number of test matrizes
nn = 100

 @testset "Test constructors standard matrix type" begin

    for iii=1:nn

        # take random dimension
        dim = rand(2:100,1,1)[1]
        # create random sparse matrix
        A = rand(dim,dim)
        A = A+A'
        # create graph from A
        g = Graph(A)

        # println(checkGraphNodes(g,A))
        @test checkGraphNodes(g,A) == true
    end
 end


 @testset "Test constructors sparse matrix type" begin

    for iii=1:nn

        # take random dimension
        dim = rand(2:100,1,1)[1]
        # create random sparse matrix
        A = sprand(dim,dim,0.1)
        A = A+A'
        # create graph from A
        g = Graph(A)

        # println(checkGraphNodes(g,A))
        @test checkGraphNodes(g,A) == true
    end
 end






