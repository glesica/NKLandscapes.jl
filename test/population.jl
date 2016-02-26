using FactCheck

import NKLandscapes
const NK = NKLandscapes
const count = 5

srand(1)

# TODO: Links gets set up properly on instantiation.

facts("population.jl") do
  context("NKLandscape") do
    l = NK.NKLandscape(4, 1)
    p = rand(NK.Population,l,count)
    @fact NK.popsize(p) --> count
    @fact length(NK.popfits(p)) --> count
  end

  context("NKqLandscape") do
    l = NK.NKqLandscape(4, 1, 2)
    p = rand(NK.Population,l,count)
    @fact NK.popsize(p) --> count
    @fact length(NK.popfits(p)) --> count
  end

  context("NKpLandscape") do
    l = NK.NKpLandscape(4, 1, 0.90)
    p = rand(NK.Population,l,count)
    @fact NK.popsize(p) --> count
    @fact length(NK.popfits(p)) --> count
  end
end
