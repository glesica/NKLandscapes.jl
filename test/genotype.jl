srand(1)

# TODO: Contribs gets updated properly when fitness(...) is called.

context("genotype.jl") do
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

    context("$(ls)") do
      context("Genotype") do
        g = NK.Genotype(0b0000, ls)
        @fact g.alleles --> 0b0000
        @fact g.landscape --> ls
      end

      context("fitness(...)") do
        g = NK.Genotype(0b0000, ls)
        @fact NK.fitness(g) --> less_than_or_equal(1.0)
        @fact NK.fitness(g) --> greater_than_or_equal(0.0)
        @fact NK.fitness(g) --> roughly(NK.contribs(g) |> mean)
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

