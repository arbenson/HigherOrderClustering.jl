# HigherOrderClustering: higher-order clustering coefficients for networks

This package provides an interface to compute "higher-order clustering coefficients" in networks, an idea developed in the following paper:

- Hao Yin, Austin R. Benson, and Jure Leskovec. [Higher-order clustering in networks](http://www.cs.cornell.edu/~arb/papers/higher-order-clustering-PRE-2018.pdf). Physical Review E, 2018.

The key computational kernel for higher-order clustering coefficients is clique enumeration. We use the [Chiba and Nishizeki algorithm](http://www.ecei.tohoku.ac.jp/alg/nishizeki/sub/j/DVD/PDF_J/J053.pdf) for this task. The package also has an interface for clique counting with this algorithm (see the examples below).

## Installation
From Julia
```julia
Pkg.add("HigherOrderClustering")
Pkg.test("HigherOrderClustering")
using HigherOrderClustering
```

## Examples
The following examples assume that you have added the package and run the following in Julia.
```julia
using HigherOrderClustering
```

#### Compute third-order clustering coefficients of C. elegans
```julia
A = load_example_data("celegans.txt")
ccfs = higher_order_ccfs(A, 3)
ccfs.global_hoccf
ccfs.avg_hoccf
ccfs.local_hoccfs
```

A "zero" value for the local clustering coefficient can mean two things. Either (i) the node is at the center of at least one higher-order wedge and none of them close or (ii) the node is not at the center of a wedge.  The field value avg_hoccf does not count the latter case in its mean.  It is easy to find the nodes in the latter case because the data structure also returns the wedge counts.
```julia
find(ccfs.ho_wedge_counts .== 0)
```

The returned data structure also includes the average clustering coefficient that considers nodes not participating in any wedges to have "0 clustering".
```julia
ccfs.avg_hoccf2
ccfs.avg_hoccf2 ≈ mean(ccfs.local_hoccfs)  # should be true
```

As the order of the clustering coefficient increases, the number of nodes that are not at the center of at least one wedge can only go up. This is something to keep in mind for analysis. When analyzing higher-order clustering, it is always useful to report the fraction of nodes that are in at least one wedge. This is about 93% for the third-order clustering coefficient in C. elegans.
```julia
sum(ccfs.ho_wedge_counts .> 0) / length(ccfs.ho_wedge_counts)
```

I recommend this paper by Marcus Kaiser that discusses the issue of nodes that do not participate in wedges (the analysis is for the classical clustering coefficient, so this set of nodes are isolated nodes and leafs):

- Marcus Kaiser. [Mean clustering coefficients: the role of isolated nodes and leafs on clustering measures for small-world networks](http://iopscience.iop.org/article/10.1088/1367-2630/10/8/083042). New Journal of Physics 10.8 (2008): 083042.



#### Clique counting in arxiv-AstroPh
```julia
A = load_example_data("arxiv-AstroPh.txt")
cliques3 = kcliques(A, 3)
cliques4 = kcliques(A, 4)
mean(cliques4), std(cliques4)
cliques5 = kcliques(A, 5)
```
Note that the reported progress will be slower for the first nodes that are processed.

#### Reading data from a text file.
There is an interface to read a graph from a text file and convert it to an adjacency matrix.

First, let's create a simple graph text file with 4 edges. In Julia:
```julia
open("testgraph.txt", "w") do f
	write(f, "# My graph file\n% is here\n1 2\n1 3\n2 3\n4 8\n")
end
```

From the command line:
```
$ cat testgraph.txt
# My graph file
% is here
1 2
1 3
2 3
4 8
```

We can read the file in one of two ways. In the first way, the indices in the matrix are kept to match the indices in the text file:
```julia
(A, inds) = read_undir_graph_txt("testgraph.txt")
size(A)  # should be (8, 8)
```
(Note that the function ignored the lines starting with '#' and '%'.)

While the adjacency matrix is 8 x 8, there are really only 5 nodes (there are no nodes labeled 5, 6, or 7). We can alternatively read the graph so that the nodes are indexed 1, 2, ..., n, where n is the total number of nodes. Setting the second argument to true does this:
```julia
(A, inds) = read_undir_graph_txt("testgraph.txt", true)
size(A)  # should be (5, 5)
```
The vector inds contains the indexing information. The entry inds[v] is the original label of the node:
```julia
inds[5]  # should be 8
```
