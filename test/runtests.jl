using HigherOrderClustering
using Base.Test
using Combinatorics

function test1()
    I = [1, 1, 1, 2, 2, 3]
    J = [2, 3, 4, 3, 4, 4]
    A = convert(SparseMatrixCSC{Int64,Int64},
                sparse(I, J, ones(length(I)), 4, 4))
    A = max.(A, A')

    @test kcliques(A, 2) == [3, 3, 3, 3]
    @test kcliques(A, 3) == [3, 3, 3, 3]
    @test kcliques(A, 4) == [1, 1, 1, 1]    
    @test kcliques(A, 5) == [0, 0, 0, 0]

    ccf2 = higher_order_ccfs(A, 2)
    @test ccf2.global_hoccf == 1.0
    @test ccf2.avg_hoccf    == 1.0
    @test ccf2.avg_hoccf2   == 1.0

    ccf3 = higher_order_ccfs(A, 3)
    @test ccf3.global_hoccf == 1.0
    @test ccf3.avg_hoccf    == 1.0
    @test ccf3.avg_hoccf2   == 1.0

    push!(I, 3, 4, 4)
    push!(J, 7, 7, 4)
    A = convert(SparseMatrixCSC{Int64,Int64},
                sparse(I, J, ones(length(I)), 7, 7))
    A = max.(A, A')

    @test sum(kcliques(A, 2)) == 8 * 2
    @test sum(kcliques(A, 3)) == 5 * 3
    @test sum(kcliques(A, 4)) == 1 * 4
    @test sum(kcliques(A, 5)) == 0 * 5
    
    @test kcliques(A, 4) == [1, 1, 1, 1, 0, 0, 0]
end

function test2()
    # Test k-cliques on a small network
    A = load_example_data("celegans.txt")
    degs = vec(sum(A, 2))
    n = size(A, 2)
    deg_order = zeros(Int64, n)
    deg_order[sortperm(vec(sum(A, 1)))] = collect(1:n)

    function higher_neighbors(j::Int64)
        oj = deg_order[j]
        nbrs = find(A[:, j])
        return filter(x -> deg_order[x] > oj, nbrs)
    end

    counts = zeros(Int64, size(A, 2))
    for j = 1:n
        for w in higher_neighbors(j)
            counts[[j, w]] += 1
        end
    end
    counts == kcliques(A, 2)
        
    counts = zeros(Int64, size(A, 2))    
    for j = 1:n
        for (w, x) in combinations(higher_neighbors(j), 2)
            if A[w, x] > 0
                counts[[j, w, x]] += 1
            end
        end
    end
    @test counts == kcliques(A, 3)

    counts = zeros(Int64, size(A, 2))
    for j = 1:n
        for (w, x, y) in combinations(higher_neighbors(j), 3)
            if A[w, x] > 0 && A[w, y] > 0 && A[x, y] > 0
                counts[[j, w, x, y]] += 1
            end
        end
    end
    @test counts == kcliques(A, 4)

    counts = zeros(Int64, size(A, 2))    
    for j = 1:n
        for (w, x, y, z) in combinations(higher_neighbors(j), 4)
            if A[w, x] > 0 && A[w, y] > 0 && A[x, y] > 0 &&
                A[w, z] > 0 && A[x, z] > 0 && A[y, z] > 0
                counts[[j, w, x, y, z]] += 1                
            end
        end
    end
    @test counts == kcliques(A, 5)
end

test1()
test2()

