using FactCheck

import NKLandscapes
const NK = NKLandscapes

const psize = 5

srand(1)

# TODO: Links gets set up properly on instantiation.

facts("population.jl") do
  landscapes = [
    NK.NKLandscape(4, 1),
    NK.NKqLandscape(4, 1, 2),
    NK.NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscapes
    pr = rand(NK.Population, ls, psize)
    pz = zeros(NK.Population, ls, psize)

    context("$(ls)") do
      context("popsize(...)") do
        @fact NK.popsize(pr) --> psize
        @fact NK.popsize(pz) --> psize
      end

      context("popfits(...)") do
        @fact length(NK.popfits(pr)) --> psize
        @fact length(NK.popfits(pz)) --> psize
      end
    end
  end
end
