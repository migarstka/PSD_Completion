module GraphModule
export Graph, numberOfVertizes,mcsSearch, mcsmSearch, findHigherNeighbors, findParent, isPerfectOrdering, isConnected


    # -------------------------------------
    # TYPE DEFINITIONS
    # -------------------------------------
    # TODO: Rename attributes in a more consistent way
    type Graph
        adjacencyList::Array{Array{Int64,1}}
        ordering::Array{Int64} # σ(v)
        reverseOrder::Array{Int64} #σ^(-1)(i)

         #constructor for adjacencylist input
        function Graph(adjacencyList::Array)
            ordering = collect(1:size(adjacencyList,1))
            g = new(adjacencyList,ordering)

            reverseOrder = zeros(size(g.ordering,1))
            # also compute reverse order σ^-1(v)
            for i = 1:size(reverseOrder,1)
                reverseOrder[Int64(g.ordering[i])] = i
            end
            g.reverseOrder = reverseOrder
            return g
        end

        # constructor for input matrix (dense data type)
        function Graph(A::Array{Float64})
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
            g = new(adjacencyList,ordering)
            # make sure that the graph is connected
            connectGraph(g)
            return g
        end

        function Graph(A::Array{Int64})
            return Graph(float(A))
        end

        # constructor for input matrix (sparse data type)
        function Graph(A::SparseMatrixCSC)
            if A != A'
                error("Please input a symmetric matrix.")
            end
            N = size(A,1)
            ordering = collect(1:1:N)
            adjacencyList = [Int64[] for i=1:N]
            for col = 1:N-1
                col_start = A.colptr[col]
                col_end = A.colptr[col+1]-1

                for i=col_start:1:col_end
                        row = A.rowval[i]
                        if row > col
                            push!(adjacencyList[row],col)
                            push!(adjacencyList[col],row)
                        end
                end
            end
            g = new(adjacencyList,ordering)
            connectGraph(g)
            return g
        end

    end #TYPE definition

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
        paths = [] #saves for each vertex the paths that were tried
        for i = N:-1:1
            # find unvisited vertex of maximum weight
            unvisited_weights = weights.*unvisited
            v = indmax(unvisited_weights)
            perfectOrdering[v] = i
            unvisited[v] = 0
            #println(" >>> Pick next vertex: v = $(v) and assign order i=$(i)\n")
            # find all unvisited vertices u with a path u, x1, x2, ..., v in G, s.t. w(xi) < w(u) and put them in set S
            # in the first step there will be no valid path, therefore choose S to be the direct neighbors
            if i == N
                S = filter(j->unvisited[j]==1,g.adjacencyList[v])
            else
                anchestors = [v]
                S = Array{Int64}[]
                S = reachableVertices(S,g,v,1,copy(unvisited),weights,anchestors)
               # println("S=$(S) of vertex $(v)")
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
        for edge in F
            push!(g.adjacencyList[edge[1]],edge[2])
            push!(g.adjacencyList[edge[2]],edge[1])
        end
        return nothing
    end

    function reachableVertices(r, g,v,depth,unvisited,weights,anchestors)
        # initialize set of reachable vertices (direct or indirect neighbors)
        # n: neighbors of v
        doPrint = false
        doPrint && println("    "^depth * "reachableVertices(v=$(v),depth=$(depth),unvisited=$(unvisited),weights=$(weights))\n")
        unvisitedNeighbors = filter(i->unvisited[i]==1,g.adjacencyList[v])
        if size(unvisitedNeighbors,1) == 0
            doPrint && println("    "^depth *"End reached at vertex $(v)\n")
            return r
        end


        # TODO: Check if necessary to sort neighbors
        # sort unvisitedNeighbors based on their weight
        sort!(unvisitedNeighbors, by=x->weights[x])
        doPrint &&  println("    "^depth *"before unvisitedNeighbors=$(unvisitedNeighbors), r=$(r)\n")

        # only visit neighbors that havent been added to r yet
        filter!(f->!in(f,r),unvisitedNeighbors)

        # direct unvisited neighbors are always added
        if depth == 1
            r = vcat(r,unvisitedNeighbors)
        end



        doPrint &&  println("    "^depth *"after unvisitedNeighbors=$(unvisitedNeighbors)\n")

        for w in unvisitedNeighbors
            doPrint &&  println("    "^depth *"Check w=$(w) of v=$(v), weights(w)=$(weights[w]),weights(v)=$(weights[v])\n")

            # check again here since might been already checked at a lower depth, i.e. unvisited variable was modified
            unvisited[w] = 0

            # case that we have a sequence w(xi) < w(w), w(target)
            if depth > 1
                doPrint &&  println("    "^depth *"validpath before w=$(w), weights(w)=$(weights[w]),anchestors=$(anchestors)\n")
                if validPath(w,weights,anchestors)
                    # prevent adding a vertex twice if he has two valid paths
                    if !in(w,r)
                        push!(r,w)
                    end
                end
            end

            # else
            #     if depth > 1 #or maybe > 2, means we tried this particular path for this vertex
            #         push!(paths[v],anchestors)
            #     end
            #end
            # sequence of vertices to reach w
            anchestorsOfW = copy(anchestors)
            push!(anchestorsOfW,w)
            doPrint && println("    "^depth *"in, r=$(r), w=$(w), unvisited=$(unvisited)\n")
            r = reachableVertices(r,g,w,depth+1,copy(unvisited),weights,anchestorsOfW)
            doPrint && println("    "^depth *"out, r=$(r)\n")
        end
        return r
    end

    function validPath(w::Int64,weights::Array{Float64},anchestors::Array{Int64})
        # trivial case of direct neighbor
        if size(anchestors,1) == 1
            return true

        # check condition max(w(xi)) < w(w),w(v) for a path w,x1,x2,...,v
        else
            v = anchestors[1]
            xs = anchestors[2:end]
            maxWeight = maximum(weights[xs])
            #println("    "^size(anchestors,1) *"v=$(v), weights(w)=$(weights[w]), xs=$(xs), maxWeight=$(maxWeight)\n")

            if weights[w] > maxWeight && weights[v] > maxWeight
                return true
            else
                return false
            end
        end

    end


    # returns lists of vertices that form the unconnected subgraphs (breath-first-search style)
    function getConnectedParts(g::Graph)
        N = numberOfVertizes(g)
        subgraphs = []
        visited = zeros(N)
        allVisited = false

         while !allVisited
            frontier = [findfirst(x->x== 0,visited)]
            visitedNodes = [frontier[1]]
            visited[frontier[1]] = 1
             while size(frontier,1) > 0
                nextFrontier = []
                for u in frontier
                    for v in g.adjacencyList[u]
                        if visited[v] == 0
                            push!(visitedNodes,v)
                            visited[v] = 1
                            push!(nextFrontier,v)
                        end
                    end
                end
                frontier = nextFrontier
            end
            # add vertices of subgraph to array
            push!(subgraphs,visitedNodes)

            # if all vertices are processed break
            if !in(0,visited)
                allVisited = true
            end
        end
        return subgraphs
    end

    # connects an unconnected graph by adding edges
    function connectGraph(g::Graph)
        subgraphs = getConnectedParts(g)

        # if more than one subgraph are found, add one edge between the first node of each subgraph
        if size(subgraphs,1) > 1
            for i=1:size(subgraphs,1)-1
                node_subgraphA = subgraphs[i][1]
                node_subgraphB = subgraphs[i+1][1]
                push!(g.adjacencyList[node_subgraphA],node_subgraphB)
                push!(g.adjacencyList[node_subgraphB],node_subgraphA)
            end
        end
        return nothing

    end

    # check if a graph is connected
    function isConnected(g::Graph)
        return size(getConnectedParts(g),1) == 1
    end
    # check if the ordering of the graph is a perfect elimination ordering (i.e. for every v, are all higher neighbors adjacent?)
    # start with lowest-order vertex v, find lowest neighbor u of v with higher order. Then verify that w is adjacent to all higher order neighbors of v
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