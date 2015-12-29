using NKLandscapes
using FactCheck

srand(1)

facts("NKLandscapes.jl") do
  landscapes = [
    NKLandscape(10, 1),
    NKLandscape(10, 2),
    NKqLandscape(10, 1, 2),
    NKpLandscape(10, 1, 0.8),
  ]

  for l = landscapes
    context("$(l)") do
      context("Landscapes") do
        g = rand(Genotype, l)
        f = fitness(g, l)


        context("Neighbors should differ at one locus") do
          function test_neighbors(genotype, landscape)
            nbrs = all_neighbors(genotype, landscape)
            for i = 1:number_neighbors(genotype, landscape)
              nbr = nbrs[:, i]
              @fact (genotype - nbr) |> sum |> abs --> 1
            end
          end

          test_neighbors(g, l)
        end

        if typeof(l) != NKLandscape
          context("Neutral neighbors should have the same fitness") do
            function test_neutral_neighbors(genotype, landscape)
              nbrs = neutral_neighbors(genotype, landscape)
              score = fitness(genotype, landscape)
              @fact size(nbrs)[1] --> landscape.n
              @fact size(nbrs)[2] --> greater_than(0)
              for j = 1:size(nbrs)[2]
                @fact fitness(nbrs[:,j], landscape) --> score
              end
            end

            test_neutral_neighbors(g, l)
          end
        end

        context("Fitter neighbors should all be fitter") do
          function test_fitter_neighbors(nbrs, landscape, score)
            @fact size(nbrs)[1] --> landscape.n
            @fact size(nbrs)[2] --> greater_than(0)
            for j = 1:size(nbrs)[2]
              @fact fitness(nbrs[:,j], landscape) --> greater_than(score)
            end
          end

          fn = fitter_neighbors(g, l)
          test_fitter_neighbors(fn, l, f)
        end

        context("Fittest n neighbors should be all neighbors") do
          fnn = fittest_neighbors(g, l, l.n)
          @fact size(fnn)[1] --> l.n
          @fact size(fnn)[2] --> number_neighbors(g, l)
        end

        context("Fittest 1 neighbor should be the fittest neighbor") do
          nbrs = fitter_neighbors(g, l)
          nbr1 = fittest_neighbor(g, l)
          @fact nbr1 --> nbrs[:,end]
        end

        context("Adaptive walks") do
          function test_adaptive_walk(walk_function, genotype, landscape)
            walk = walk_function(genotype, landscape)
            @fact walk.length --> greater_than(0)
            for i = 2:walk.length
              @fact fitness(walk.history[:,i], landscape) --> greater_than(fitness(genotype, landscape))
            end
          end

          context("Random adaptive walk should terminate and move uphill") do
            test_adaptive_walk(random_walk, g, l)
          end

          context("Greedy adaptive walk should terminate and move uphill") do
            test_adaptive_walk(greedy_walk, g, l)
          end

          context("Reluctant adaptive walk should terminate and move uphill") do
            test_adaptive_walk(reluctant_walk, g, l)
          end
        end
      end

      context("Populations") do
        context("Should be the correct size") do
          function test_population_size(p::Population, l::Landscape, n)
            @fact popsize(p) --> n
            @fact size(p)[1] --> l.n
            @fact size(p)[2] --> n
          end

          for n = [1, 10, 100]
            rp = rand(Population, l, n)
            test_population_size(rp, l, n)
            zp = zeros(Population, l, n)
            test_population_size(zp, l, n)
          end
        end

        context("Should compute fitnesses") do
          function test_population_fitnesses(p::Population, l::Landscape)
            fs = popfits(p, l)
            for i = 1:popsize(p)
              @fact fs[i] --> fitness(p[:,i], l)
            end
          end

          for n = [1, 10, 100]
            rp = rand(Population, l, n)
            test_population_fitnesses(rp, l)
            zp = zeros(Population, l, n)
            test_population_fitnesses(zp, l)
          end
        end
      end

      for n = [1, 10, 100]
        context("Selection") do
          context("Should result in a population of the same size") do
            rp = rand(Population, l, n)
            np = propsel(rp, l)
            @fact popsize(np) --> n
          end

          context("Should contain only members of the original population") do
            rp = rand(Population, l, n)
            np = propsel(rp, l)
            for i_n = 1:n
              founds = map(1:n) do i_o
                np[:,i_n] == rp[:,i_o]
              end
              @fact any(founds) --> true
            end
          end

          context("Should create a new population when using propsel") do
            rp = rand(Population, l, n)
            np = propsel(rp, l)
            @fact is(rp, np) --> false
          end
        end

        context("Mutation") do
          context("Should result in a population of the same size") do
            rp = rand(Population, l, n)
            np = bwmutate(rp, l, 1.0)
            @fact popsize(np) --> n
          end

          context("Should mutate population members") do
            rp = rand(Population, l, n)
            np = bwmutate(rp, l, 1.0)
            @fact np --> not(rp)
          end
        end
      end
    end
  end
end
