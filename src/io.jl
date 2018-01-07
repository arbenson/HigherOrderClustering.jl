"""
Read a graph from a text file. Skips any lines that begin with '#' or '%'
characters. If oneindex is true, then nodes are mapped to values 1, ..., n.

Returns (A, v), where
   A is the n x n undirected graph adjacency matrix
   v is a length-n vector such that v[i] is the original node identifier
     (if oneindex is false, then v[i] = i, i = 1, ..., n)
"""
function read_undir_graph_txt(filename::AbstractString, oneindex::Bool=false)
    # index mapping
    index_map = Dict{Int64,Int64}()
    index_map_vec = Int64[]
    function get_mapped_index(x::Int64)
        if !haskey(index_map, x)
            next_index = length(index_map) + 1
            index_map[x] = next_index
            push!(index_map_vec, x)
            return next_index
        end
        return index_map[x]
    end

    # Read data
    I = Int64[]
    J = Int64[]
    open(filename) do f
        for line in eachline(f)
            # Skip lines starting with '#' or '%'
            if line[1] == '#' || line[1] == '%'; continue; end
            edge = split(line)
            u = parse(Int64, edge[1])
            v = parse(Int64, edge[2])
            if oneindex
                push!(I, get_mapped_index(u))
                push!(J, get_mapped_index(v))
            else
                push!(I, u)
                push!(J, v)
            end
        end
    end

    # Form adjacency matrix
    if min(minimum(I), minimum(J)) < 1
        error("Minimum node value is less than 1. Try setting ondeindex parameter to true.")
    end
    n = max(maximum(I), maximum(J))
    A = convert(SparseMatrixCSC{Int64,Int64}, sparse(I, J, ones(length(I)), n, n))
    A = max.(A, A')

    if oneindex; return (A, index_map_vec); end
    return (A, collect(1:n))
end

"""
Load an example file from the data directory.
"""
function load_example_data(filename::AbstractString)
    pathname = joinpath(dirname(dirname(@__FILE__)), "data")
    filename = joinpath(pathname, filename)
    if   isfile(filename); return read_undir_graph_txt(filename, true)[1]
    else error(@sprintf("Could not find file %s", name))
    end
end
