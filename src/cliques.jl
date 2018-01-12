using MatrixNetworks

function _kcliques(A::SparseMatrixCSC{Int64,Int64}, k::Int64)
    A = min.(A, 1)
    A -= spdiagm(diag(A))    
    # "global" variables for the recursive algorithm
    n = size(A, 2)
    adjacency_lists = Vector{Vector{Int64}}(n)
    for v = 1:n
        adjacency_lists[v] = filter(x -> x != v, find(A[:, v]))
    end
    C = Int64[]
    labels = k * ones(Int64, n)
    counts = zeros(Int64, n)

    # Get all neighbors of a node with a given label
    function neighbors(node::Int64, label::Int64)
        index = 1
        for nbr in adjacency_lists[node]
            nbr_label = labels[nbr]
            if nbr_label == label; index += 1; end
            if nbr_label  > label; break; end
        end
        return adjacency_lists[node][1:(index - 1)]
    end

    # Get the number of neighbors of a node with a given label
    function num_neighbors(node::Int64, label::Int64)
        index = 1
        for nbr in adjacency_lists[node]
            nbr_label = labels[nbr]
            if nbr_label == label; index += 1; end
            if nbr_label  > label; break; end
        end
        return index - 1
    end

    # Get the nodes in U in sorted order that need to be processed
    function nodes_to_process(U::Vector{Int64}, label::Int64)
        sorted_nbrs = sort([(num_neighbors(v, label), v) for v in U], rev=true)
        return [x[2] for x in sorted_nbrs if x[1] > 1]
    end

    # Adjust adjacency lists of neighbors of nodes in U to put neighbors with a
    # given label in the front.
    function move_front(U::Vector{Int64}, label::Int64)
        for v in U
            front_nbrs = Int64[]
            back_nbrs  = Int64[]
            for nbr in adjacency_lists[v]
                if labels[nbr] == label; push!(front_nbrs, nbr)
                else                     push!(back_nbrs, nbr)
                end
            end
            num_front = length(front_nbrs)
            total = length(adjacency_lists[v])
            adjacency_lists[v][1:num_front] = front_nbrs
            adjacency_lists[v][(num_front + 1):total] = back_nbrs
        end
    end

    # Adjust adjacency lists of neighbors of nodes in U to put node vi directly
    # after nodes with a given label.
    function move_after(U::Vector{Int64}, vi::Int64, label::Int64)
        for v in U
            front_nbrs = Int64[]
            for nbr in adjacency_lists[v]
                if nbr == vi; continue; end
                if labels[nbr] <= label; push!(front_nbrs, nbr)
                else                     break
                end
            end
            num_front = length(front_nbrs)
            adjacency_lists[v][1:num_front] = front_nbrs
            adjacency_lists[v][num_front + 1] = vi
        end
    end

    # Process all current cliques, which each consist of:
    #   1. the nodes in the global stack C
    #   2. a node in U
    #   3. a neighbor of u with label 2
    function process_cliques(U::Vector{Int64})
        C_cnt = 0
        for v in U
            v_cnt = 0
            for nbr in neighbors(v, 2)
                # Only count (nbr, v) edge once.
                if nbr > v
                    counts[nbr] += 1
                    v_cnt += 1
                end
            end
            counts[v] += v_cnt
            C_cnt += v_cnt
        end
        counts[C] += C_cnt
    end

    # Recursive algorithm of Chiba and Nishizeki.
    function rcliques(U::Vector{Int64}, r::Int64)
        # Check base case
        if r == 2
            process_cliques(U)
            return
        end
        
        for (ind, v) in enumerate(nodes_to_process(U, r))
            if r == k
                print(@sprintf("%d of %d (%d-cliques) \r", ind, length(U), r))
                flush(STDOUT)
                if ind == length(U); println(""); end
            end
            # Get neighborhood of v
            Up = neighbors(v, r)
            labels[Up] = r - 1
            move_front(Up, r - 1)

            # Recurse on neighborhood of node v
            push!(C, v)
            rcliques(Up, r - 1)
            pop!(C)

            # Restore neighborhood
            labels[Up] = r

            # Eliminate v
            labels[v] = r + 1
            move_after(Up, v, r)
        end
    end

    rcliques(find(vec(sum(A, 2)) .>= (k - 1)), k)
    return counts
end

"""
Count k-cliques in the graph A using the Chiba-Nishizeki algorithm.

Returns a vector c, where c[v] is the number of k-cliques containing node v
in the graph.
"""
function kcliques(A::SparseMatrixCSC{Int64,Int64}, k::Int64)
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    clique_counts = zeros(Int64, size(A, 1))
    (d, rt) = corenums(A)
    inds = find(d .>= (k - 1))
    counts = _kcliques(A[inds, inds], k)
    clique_counts[inds] = counts
    return clique_counts
end
