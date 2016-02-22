using FactCheck

import NKLandscapes
const NK = NKLandscapes

srand(1)

facts("landscape.jl") do
  context("NKLandscape") do
    l = NK.NKLandscape(4, 1)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n
  end

  context("NKqLandscape") do
    l = NK.NKqLandscape(4, 1, 2)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n
  end

  context("NKpLandscape") do
    l = NK.NKpLandscape(4, 1, 0.90)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n
  end
end
