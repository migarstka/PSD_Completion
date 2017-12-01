module GraphModule
export Graph, numberOfVertizes,mcsSearch, mcsmSearch, findHigherNeighbors, findParent, isPerfectOrdering


    # -------------------------------------
    # TYPE DEFINITIONS
    # -------------------------------------
    # TODO: Rename attributes in a more consistent way
    type Graph
        adjacencyList::Array{Array{Int64,1}}
        ordering::Array{Int64}
        reverseOrder::Array{Int64}

         #constructor for adjacencylist input
        function Graph(adjacencyList::Array{Array{Float64}})
            ordering = collect(1:size(adjacencyList,1))
            new(adjacencyList,ordering)
        end

        # constructor for input matrix (dense data type)
        function Graph(A::Array{Float64})
            println("Es")
            if A != A'
                error("Please input a symmetric matrix.")
            end
            N = size(A,1)
            ordering = collect(1:1:N)
            adjacencyList = [Int64[] for i=1:N]
            for j = 1:N-1
                for i=j+1:N
                    if A[i,j] != 0
                        push!(adjacencyList[i],j)
                        push!(adjacencyList[j],i)
                    end
                end
            end
            new(adjacencyList,ordering)
        end

        # constructor for input matrix (sparse data type)
        function Graph(A::Array{Float64})
            println("Es")
            if A != A'
                error("Please input a symmetric matrix.")
            end
            N = size(A,1)
            ordering = collect(1:1:N)
            adjacencyList = [Int64[] for i=1:N]
            for j = 1:N-1
                for i=j+1:N
                    if A[i,j] != 0
                        push!(adjacencyList[i],j)
                        push!(adjacencyList[j],i)
                    end
                end
            end
            new(adjacencyList,ordering)
        end
    end

    # -------------------------------------
    # FUNCTION DEFINITIONS
    # -------------------------------------
    function numberOfVertizes(g::Graph)
        return size(g.ordering,1)
    end

    # returns the neighbor with the lowest order higher than the nodes order
    function findParent(g::Graph,higherNeighbors::Array{Int64})
        if size(higherNeighbors,1) > 0
            return higherNeighbors[indmin(g.ordering[higherNeighbors])]
        else
            return 0
        end
    end

    function findHigherNeighbors(g::Graph,nodeNumber::Int64)
        order = g.ordering[nodeNumber]
        neighbors = g.adjacencyList[nodeNumber]
        higherNeighbors = neighbors[find(f->f>order,g.ordering[neighbors])]
        return higherNeighbors
    end

        function findLowerNeighbors(g::Graph,nodeNumber::Int64)
        order = g.ordering[nodeNumber]
        neighbors = g.adjacencyList[nodeNumber]
        lowerNeighbors = neighbors[find(f->f<order,g.ordering[neighbors])]
        return lowerNeighbors
    end

    # performs a maximum cardinality search and updates the ordering to the graph (only perfect elim. ordering if graph is chordal)
    function mcsSearch(g::Graph)
        N = numberOfVertizes(g)
        weights = zeros(N)
        unvisited = ones(N)
        perfectOrdering = zeros(N)
        for i = N:-1:1
            # find unvisited vertex of maximum weight
            unvisited_weights = weights.*unvisited
            indMax = indmax(unvisited_weights)
            perfectOrdering[indMax] = i
            unvisited[indMax] = 0
            for neighbor in g.adjacencyList[indMax]
                if unvisited[neighbor] == 1
                    weights[neighbor]+=1
                end
            end
        end
        # update ordering of graph
        g.ordering = perfectOrdering
        reverseOrder = zeros(size(perfectOrdering,1))
        # also compute reverse order σ^-1(v)
        for i = 1:N
            reverseOrder[Int64(perfectOrdering[i])] = i
        end
        g.reverseOrder = reverseOrder
        return nothing
    end

    # implementation of the MCS-M algorithm (see. A. Berry - Maximum Cardinality Search for Computing Minimal Triangulations of Graphs) that finds a minimal triangulation
    function mcsmSearch(g::Graph)
        # initialize edge set F of fill-in edges
        F = []
        N = numberOfVertizes(g)
        weights = zeros(N)
        unvisited = ones(N)
        perfectOrdering = zeros(N)
        for i = N:-1:1
            # find unvisited vertex of maximum weight
            unvisited_weights = weights.*unvisited
            v = indmax(unvisited_weights)
            perfectOrdering[v] = i
            unvisited[v] = 0

            # find all unvisited vertices u with a path u, x1, x2, ..., v in G, s.t. w(xi) < w(u) and put them in set S
            # in the first step there will be no valid path, therefore choose S to be the direct neighbors
            if i == N
                S = filter(i->unvisited[i]==1,g.adjacencyList[v])
            else
                S = reachableVertices(g,v,0,1,copy(unvisited),weights)
            end
            # increment weight of all vertices w and if w and v are no direct neighbors, add edges to F
            for w in S
                weights[w]+=1
                if !in(w,g.adjacencyList[v])
                    push!(F,[w,v])
                end
            end
        end

        # update ordering of graph
        g.ordering = perfectOrdering
        # also compute reverse order σ^-1(v)
        reverseOrder = zeros(size(perfectOrdering,1))
        for i = 1:N
            reverseOrder[Int64(perfectOrdering[i])] = i
        end
        g.reverseOrder = reverseOrder

        # TODO: Do other algorithms break if the adjacencylist is not ordered anymore? -> if so: sort adjacencylist
        # make graph chordal by adding the new edges of F to E
        @show(F)
        for edge in F
            push!(g.adjacencyList[edge[1]],edge[2])
            push!(g.adjacencyList[edge[2]],edge[1])
        end
        return nothing
    end

    function reachableVertices(g,v,refweight,depth,unvisited,weights)
        # initialize set of reachable vertices (direct or indirect neighbors)
        r = []
        # n: neighbors of v
        unvisitedNeighbors = filter(i->unvisited[i]==1,g.adjacencyList[v])
        if size(unvisitedNeighbors,1) == 0
            return r
        end

        for w in unvisitedNeighbors
            # check again here since might been already checked at a lower depth, i.e. unvisited variable was modified
            if unvisited[w] == 1
                unvisited[w] = 0
                if depth == 1
                    r = vcat(r,w)
                    refweight = weights[w]
                    rv= reachableVertices(g,w,refweight,depth+1,unvisited,weights)
                else
                    if weights[w] <= weights[v]
                        rv= reachableVertices(g,w,refweight,depth+1,unvisited,weights)
                    else #w(w) > w(v)
                        if weights[w] > refweight
                            r = vcat(r,w)
                            rv = reachableVertices(g,w,refweight,depth+1,unvisited,weights)
                        else
                          refweight = weights[w]
                          rv=  reachableVertices(g,w,refweight,depth+1,unvisited,weights)
                        end
                    end
                end
                r = vcat(r,rv)
            else
                if depth == 1
                    r = union(r,w)
                end
            end
        end
        return r
    end


    # check if the ordering of the graph is a perfect elimination ordering (i.e. for every v, are all higher neighbors adjacent?)
    # start with lowest-order vertex v, find lowest neighbor w of v with higher order. Then verify that w is adjacent to all higher order neighbors of v
    # Algorithm has running time O(m+n)
    function isPerfectOrdering(g::Graph)
        for v in g.reverseOrder
            higherNeighbors = findHigherNeighbors(g,v)
            if size(higherNeighbors,1) > 1
                u = higherNeighbors[indmin(g.ordering[higherNeighbors])]
                for w in higherNeighbors
                    if w != u && !any(x->x==u,g.adjacencyList[w])
                            return false
                    end
                end
            end
        end
        return true
    end
end #MODULE