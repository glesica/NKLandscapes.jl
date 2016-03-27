srand(1)

context("walks.jl") do
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

    # Get the min fitness genotype so we know we can improve.
    gmin = NK.Genotype(0b000, ls)
    for a::NK.AlleleString = 0b000:0b1111
      g = NK.Genotype(a, ls)
      if NK.fitness(g) < NK.fitness(gmin)
        gmin = g
      end
    end

    context("$(ls)") do
      context("Adaptive walks") do
        function test_adaptive_walk(walk_function, genotype)
          walk = walk_function(genotype)
          @fact walk.length --> greater_than_or_equal(0)
          for i = 2:walk.length
            @fact NK.fitness(walk.history_list[i]) --> greater_than(NK.fitness(genotype))
          end
        end

        context("Random adaptive walk should terminate and move uphill") do
          test_adaptive_walk(NK.random_adaptive_walk, gmin)
        end

        context("Greedy adaptive walk should terminate and move uphill") do
          test_adaptive_walk(NK.greedy_adaptive_walk, gmin)
        end

        context("Reluctant adaptive walk should terminate and move uphill") do
          test_adaptive_walk(NK.reluctant_adaptive_walk, gmin)
        end
      end

      context("Neutral walks") do
        function test_neutral_walk(walk_function, genotype)
          walk = walk_function(genotype)
          @fact walk.length --> greater_than_or_equal(0)
          f = NK.fitness(genotype)
          for i = 1:walk.length
            f0 = NK.fitness(walk.history_list[i])
            @fact f0 --> roughly(f)
          end
        end

        context("Neutral walk should terminate and move sideways") do
          test_neutral_walk(NK.neutral_walk, gmin)
        end
      end
    end
  end
end

