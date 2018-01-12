module HigherOrderClustering

export hoccf_data, higher_order_ccfs, kcliques
export read_undir_graph_txt, load_example_data

include("io.jl")
include("cliques.jl")

"""
hoccf_data
----------

This is the data type returned by the function that computes higher-order clustering coefficients.

The field values are

order::Int64
    The order of the clustering coefficient (order = 2 is the classical
    clustering coefficient definition).

global_hoccf::Float64
    The global higher-order clustering coefficient.

avg_hoccf::Float64
    The average higher-order clustering coefficient (the mean is taken over
    nodes at the center of at least one wedge).

avg_hoccf2::Float64
     The average higher-order clustering coefficient, where the local clustering
    of a node not in any wedge is considered to be 0.

local_hoccfs::Vector{Float64}
    The vector of local higher-order clustering coefficients. If a node is not
    the center of at least one wedge, then the value is 0.

ho_wedge_counts::Vector{Int64}
    Higher-order wedge counts of nodes: ho_wedge_counts[v] is the number of higher-order
    wedges with node v at its center.

clique_counts::Vector{Int64}
    Clique counts of nodes: clique_counts[v] is the number of k-cliques containing
    node v where k = order + 1.
"""
immutable hoccf_data
    order::Int64
    global_hoccf::Float64
    avg_hoccf::Float64
    avg_hoccf2::Float64
    local_hoccfs::Vector{Float64}
    ho_wedge_counts::Vector{Int64}
    clique_counts::Vector{Int64}
end

"""
higher_order_ccfs
-----------------

This is the function that computes higher-order clustering coefficients.
It takes two inputs:

A::SparseMatrixCSC{Int64,Int64}
    The adjacency matrix of an undirected graph.

l::Int64
    The order of the clustering coefficient (l = 2 is the classical definition).
"""
function higher_order_ccfs(A::SparseMatrixCSC{Int64,Int64}, l::Int64)
    A = min.(A, 1)
    A -= spdiagm(diag(A))
    n = size(A, 1)
    # Get clique counts
    clique_counts1 = kcliques(A, l)
    clique_counts2 = kcliques(A, l + 1)
    degs = vec(sum(A, 2))
    # normalize
    wedge_counts = (degs - l + 1) .* clique_counts1
    nz_inds = find(wedge_counts .> 0)
    local_hoccfs = zeros(Float64, n)
    local_hoccfs[nz_inds] = l * clique_counts2[nz_inds] ./ wedge_counts[nz_inds]
    return hoccf_data(l, l * sum(clique_counts2) / sum(wedge_counts),
                      mean(local_hoccfs[nz_inds]), mean(local_hoccfs),
                      local_hoccfs, wedge_counts, clique_counts2)
end

end # module
