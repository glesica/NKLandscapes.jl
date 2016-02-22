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

    context("$(ls)") do
      context("number_neighbors(...)") do
        @fact NK.number_neighbors(g) --> 4
      end

      context("all_neighbors(...)") do
        nbrs = NK.all_neighbors(g)
      end

      context("random_neighbor(...)") do
      end

      context("neutral_neighbors(...)") do
      end

      context("fitter_neighbors(...)") do
      end

      context("fitness_range_neighbors(...)") do
      end

      context("fittest_neighbors(...)") do
      end

      context("fittest_neighbor(...)") do
      end
    end
  end
end
