using FactCheck

import NKLandscapes
const NK = NKLandscapes
const count = 5

srand(1)

# TODO: Links gets set up properly on instantiation.

facts("population.jl") do
  context("NKLandscape") do
    l = NK.NKLandscape(4, 1)
    pr = rand(NK.Population,l,count)
    pz = zeros(NK.Population,l,count)
    @fact NK.popsize(pr) --> count
    @fact NK.popsize(pz) --> count
    @fact length(NK.popfits(pr)) --> count
    @fact length(NK.popfits(pz)) --> count
  end

  context("NKqLandscape") do
    l = NK.NKqLandscape(4, 1, 2)
    pr = rand(NK.Population,l,count)
    pz = zeros(NK.Population,l,count)
    @fact NK.popsize(pr) --> count
    @fact NK.popsize(pz) --> count
    @fact length(NK.popfits(pr)) --> count
    @fact length(NK.popfits(pz)) --> count
  end

  context("NKpLandscape") do
    l = NK.NKpLandscape(4, 1, 0.90)
    pr = rand(NK.Population,l,count)
    pz = zeros(NK.Population,l,count)
    @fact NK.popsize(pr) --> count
    @fact NK.popsize(pz) --> count
    @fact length(NK.popfits(pr)) --> count
    @fact length(NK.popfits(pz)) --> count
  end
end
