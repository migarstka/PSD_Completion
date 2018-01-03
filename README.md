# Positive Semidefinite Completion Routine

In a positive semideﬁnite completion problem, one wants to fill in the entries of a partially specified dense matrix A
The diagonal entries A[ii] and the off-diagonal entries A[ij] with {i,j} ∈ E are ﬁxed, the other entries are free. The provided routine can find
values for the free entries that make the matrix A positive semidefinite.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

The code is written for Julia v0.6.


### Installing / Using

- Clone the repository to your local machine.
- Include the completion, tree and graph module into your project.
- Use psdCompletion() to complete your matrix A

```
workspace()
include("graph.jl")
include("tree.jl")
include("completion.jl")
include("helper.jl")

using Completion, Helper, Base.Test

# generate completable test matrix
rng = MersenneTwister(123554);
dim = 6
density = 0.1
A = generateCompleatableMatrix(dim,density,rng)

# perform positive definite completion with own implementation
W = psdCompletion(A)

# check if completed matrix is positive definite
@test isposdef(W) == true
```

## Authors

* **Michael Garstka**


## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The code is based on the following work
The implementation is based on the supernodal-multifrontal algorithms
described in these papers:
      
      L. Vandenberghe and M. S. Andersen, `Chordal Graphs and 
      Semidefinite Optimization <http://seas.ucla.edu/~vandenbe/publications/chordalsdp.pdf>`_, *Foundations and Trends in Optimization*, 2015. [`doi <http://dx.doi.org/10.1561/2400000006>`__ | `bib <http://www.doi2bib.org/#/doi/10.1561/2400000006>`__ ]
      
      M. S. Andersen, J. Dahl, and L. Vandenberghe, `Logarithmic barriers
      for sparse matrix cones <http://arxiv.org/abs/1203.2742>`_, *Optimization Methods and Software*, 2013. [`doi <http://dx.doi.org/10.1080/10556788.2012.684353>`__ | `bib <http://www.doi2bib.org/#/doi/10.1080/10556788.2012.684353>`__ ]
