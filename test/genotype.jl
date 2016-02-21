using FactCheck

import NKLandscapes as NK

srand(1)

facts("genotype.jl") do
  landscapes = [
    NKLandscape(4, 1),
    NKqLandscape(4, 1, 2),
    NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscape
    # The ith locus is linked to the locus i places to its right, assuming a
    # periodic boundary.
    ls.links = [
      0b1000,
      0b0100,
      0b0010,
      0b0001
    ]

    context("$(ls)") do
      context("Genotype") do
        g = Genotype(0b0000, ls)
        @fact g.alleles --> 0b0000
        @fact g.landscape --> ls
      end

      context("fitness(...)") do
        g = Genotype(0b0000, ls)
        @fact fitness(g) --> less_than_or_equal(1.0)
        @fact fitness(g) --> greater_than_or_equal(0.0)
        @fact fitness(g) --> roughly(contribs(g)[1])
      end

      context("rand(::Type{Genotype}, ...)") do
        g = rand(NK.Genotype, ls)
        @fact g.alleles --> less_than_or_equal(0b1111)
        @fact g.alleles --> greater_than_or_equal(0b0000)
      end

      context("zeros(::Type{Genotype}, ...)") do
        g = zeros(NK.Genotype, ls)
        @fact g.alleles --> 0
        @fact g.landscape --> ls
      end
    end
  end
end

