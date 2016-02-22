using FactCheck

import NKLandscapes
const NK = NKLandscapes

srand(1)

facts("neighbors.jl") do
  landscapes = [
    NK.NKLandscape(4, 1),
    NK.NKqLandscape(4, 1, 2),
    NK.NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscapes
    # The ith locus is linked to the locus i places to its right, assuming a
    # periodic boundary.
    ls.links = [
      0b1000,
      0b0100,
      0b0010,
      0b0001
    ]

    g = NK.Genotype(0b0000, ls)
    actual = [0b0001, 0b0010, 0b0100, 0b1000]

    context("$(ls)") do
      context("number_neighbors(...)") do
        @fact NK.number_neighbors(g) --> 4
      end

      context("all_neighbors(...)") do
        nbrs = NK.all_neighbors(g)
        @fact length(nbrs) --> 4
        @fact nbrs[1].alleles --> anyof(actual...)
        @fact nbrs[1].landscape --> ls
        @fact nbrs[2].alleles --> anyof(actual...)
        @fact nbrs[2].landscape --> ls
        @fact nbrs[3].alleles --> anyof(actual...)
        @fact nbrs[3].landscape --> ls
        @fact nbrs[4].alleles --> anyof(actual...)
        @fact nbrs[4].landscape --> ls
      end

      context("random_neighbor(...)") do
        rnbr = NK.random_neighbor(g)
        @fact rnbr.alleles --> anyof(actual...)
        @fact rnbr.landscape --> ls
      end

      context("neutral_neighbors(...)") do
      end

      context("fitter_neighbors(...)") do
      end

      context("fitness_range_neighbors(...)") do
      end

      context("fittest_neighbors(...)") do
        nbrs = NK.fitter_neighbors(g, 2)
        @fact length(nbrs) --> 2
        @fact nbrs[1].alleles --> anyof(actual...)
        @fact nbrs[1].landscape --> ls
        @fact nbrs[2].alleles --> anyof(actual...)
        @fact nbrs[2].landscape --> ls
      end

      context("fittest_neighbor(...)") do
        fnbr = NK.fittest_neighbor(g)
        @fact fnbr.alleles --> anyof(actual...)
        @fact fnbr.landscape --> ls
      end
    end
  end
end
