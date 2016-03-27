srand(1)

context("neighbors.jl") do
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
        # TODO: Fake out `contribs` to guarantee neighbors.
        nbrs = NK.neutral_neighbors(g)
        for nbr = nbrs
          @fact NK.fitness(nbr) --> NK.fitness(g)
          @fact nbr.alleles --> anyof(actual...)
          @fact nbr.landscape --> ls
        end
      end

      context("fitter_neighbors(...)") do
        # TODO: Fake out `contribs` to guarantee neighbors.
        nbrs = NK.fitter_neighbors(g)
        for nbr = nbrs
          @fact NK.fitness(nbr) --> greater_than(fitness(g))
          @fact nbr.alleles --> anyof(actual...)
          @fact nbr.landscape --> ls
        end
      end

      function test_fitness_range_neighbors(lb::Float64, ub::Float64)
        nbrs = NK.fitness_range_neighbors(g, lb, ub)
        if lb == 0.0 && ub > 1.0
          @fact length(nbrs) --> 4
        end
        for nbr = nbrs
          @fact NK.fitness(nbr) --> greater_than_or_equal(lb)
          @fact NK.fitness(nbr) --> less_than(ub)
          @fact nbr.alleles --> anyof(actual...)
          @fact nbr.landscape --> ls
        end
      end

      context("fitness_range_neighbors(...)") do
        test_fitness_range_neighbors(0.0, 1.0 + eps(Float64))
        test_fitness_range_neighbors(0.25, 0.75)
      end

      context("fittest_neighbors(...)") do
        nbrs = NK.fittest_neighbors(g, 2)
        allfits = [NK.fitness(x) for x = NK.all_neighbors(g)] |> sort
        @fact length(nbrs) --> 2
        for nbr = nbrs
          @fact NK.fitness(nbr) --> anyof(allfits[(end - 1):end]...)
          @fact nbr.alleles --> anyof(actual...)
          @fact nbr.landscape --> ls
        end
      end

      context("fittest_neighbor(...)") do
        fnbr = NK.fittest_neighbor(g)
        allfits = [NK.fitness(x) for x = NK.all_neighbors(g)]
        @fact NK.fitness(fnbr) --> maximum(allfits)
        @fact fnbr.alleles --> anyof(actual...)
        @fact fnbr.landscape --> ls
      end
    end
  end
end
