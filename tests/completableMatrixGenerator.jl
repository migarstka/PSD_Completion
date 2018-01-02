# Test File to check if the generated completable matrices fulfill the requirements:
# - symmetric
# - chordal sparsity pattern
# - nonzero diagonal

workspace()
include("../graph.jl")
include("../tree.jl")
include("../helper.jl")
include("../chompack_wrapper.jl")

using Helper, Base.Test, GraphModule, TreeModule,chomWrap

rng = MersenneTwister(222);



# Define number of test matrices
nn = 1000

@testset "Generate Pos-Def-Completable Matrices" begin
    for iii=1:nn
        if mod(iii,100) == 0
            println("iii=$(iii)")
        end
        dim = rand(rng,2:100)
        density = rand(rng,0.1:0.1:0.4)
        X,gOld = generateCompleatableMatrix(dim,density,rng)

        # check nonzero diagonals
        numDiagZeros = size(find(i->X[i,i] == 0,1:dim),1)
        @test numDiagZeros == 0

        #check symmetry
        @test X == X'

        #check chordality of the sparsity pattern
        g = Graph(X)
        mcsSearch!(g)
        @test isPerfectOrdering(g) == true

        # check if all clique-related submatrices are positive definite
        t = createTreeFromGraph(g)
        snet = createSupernodeEliminationTree(t,g)
        ct = createCliqueTree(snet,g)
        allCliquesPosDef = true
        for clique in ct.nodes
            indizes = clique.value_btm
            indizes = union(indizes,clique.value_top)
            if !isposdef(X[indizes,indizes])
                allCliquesPosDef = false
            end
        end
        @test allCliquesPosDef = true

        # double check with chompack cliques
        W,symb, Z = completeChompack(X)
        cliques = symb[:cliques]()
        allCliquesPosDef2 = true
        for clique in cliques
            indizes = clique + 1
            if !isposdef(X[indizes,indizes])
                allCliquesPosDef2 = false
            end
        end
        @test allCliquesPosDef2 = true

    end

end



