# HigherOrderClustering: higher-order clustering coefficients for networks

[![Build Status](https://travis-ci.org/arbenson/HigherOrderClustering.jl.svg?branch=master)](https://travis-ci.org/arbenson/HigherOrderClustering.jl)

[![Coverage Status](https://coveralls.io/repos/arbenson/HigherOrderClustering.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/arbenson/HigherOrderClustering.jl?branch=master)

[![codecov.io](http://codecov.io/github/arbenson/HigherOrderClustering.jl/coverage.svg?branch=master)](http://codecov.io/github/arbenson/HigherOrderClustering.jl?branch=master)

This package provides an interface to compute "higher-order clustering
coefficients" in networks, an idea developed in the paper

- Hao Yin, Austin R. Benson, and Jure Leskovec. Higher-order clustering in networks. [arXiv:1704.03913](https://arxiv.org/abs/1704.03913), 2017.

The key computational kernel for
higher-order clustering coefficients is clique enumeration. We use the
[Chiba and Nishizeki algorithm](http://www.ecei.tohoku.ac.jp/alg/nishizeki/sub/j/DVD/PDF_J/J053.pdf)
for this task. The package also has an interface for clique counting with this algorithm (see the examples below).

## Installation
From Julia
```
Pkg.add("HigherOrderClustering")
Pkg.test("HigherOrderClustering")
using HigherOrderClustering
```

## Examples
The following examples assume that you have added the package and run the following in Julia.
```
using HigherOrderClustering
```

#### Compute third-order clustering coefficients of C. elegans
```
A = load_example_data("celegans.txt")
ccfs = higher_order_ccfs(A, 3)
ccfs.global_hoccf
ccfs.avg_hoccf
ccfs.local_hoccfs
```

#### Clique counting in arxiv-AstroPh
```
A = load_example_data("arxiv-AstroPh.txt")
cliques3 = kcliques(A, 3)
cliques4 = kcliques(A, 4)
mean(cliques4), std(cliques4)
cliques5 = kcliques(A, 5)
```
Note that the reported progress will be slower for the first nodes that are
processed.

#### Reading data from a text file.
There is an interface to read a graph from a text file and convert it to
an adjacency matrix.

First, let's create a simple graph text file with 4 edges. In Julia:
```
open("testgraph.txt", "w") do f; write(f, "# My graph file\n% is here\n1 2\n1 3\n2 3\n4 8\n"); end;
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

We can read the file in one of two ways. In the first way, the indices in the
matrix are kept to match the indices in the text file:
```
(A, inds) = read_undir_graph_txt("testgraph.txt")
size(A)  # should be (8, 8)
```
(Note that the function ignored the lines starting with '#' and '%'.)

While the adjacency matrix is 8 x 8, there are really only 5 nodes (there are no
nodes labeled 5, 6, or 7). We can alternatively read the graph so that the nodes
are indexed 1, 2, ..., n, where n is the total number of nodes. Setting the second
argument to true does this:
```
(A, inds) = read_undir_graph_txt("testgraph.txt", true)
size(A)  # should be (5, 5)
```
The vector inds contains the indexing information. The entry
inds[v] is the original label of the node:
```
inds[5]  # should be 8
```
