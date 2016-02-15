using NKLandscapes
using FactCheck

srand(1)

facts("NKLandscapes.jl") do
  landscapes = [
    NKLandscape(10, 0),
    NKLandscape(10, 1),
    NKLandscape(10, 2),
    NKLandscape(10, 3),
    NKLandscape(10, 4),
    NKLandscape(10, 5),
    NKLandscape(10, 6),
    NKLandscape(10, 7),
    NKLandscape(10, 8),
    NKLandscape(10, 9),
    NKLandscape(10, 4; near=true),
    NKqLandscape(10, 1, 2),
    NKpLandscape(10, 1, 0.90),
  ]

  for l = landscapes
    context("$(l)") do
      context("Enumeration") do
        g = rand(Genotype, l)
        f = fitness(g, l)

        context("gtoi and itog") do
          ig = gtoi(g, l)
          gg = itog(ig, l)
          @fact gg --> g
        end

        context("fitness") do
          ig = gtoi(g, l)
          @fact fitness(ig, l) --> fitness(g, l)
        end
      end

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

        context("Random neighbor should be a neighbor") do
          g0 = random_neighbor(g, l)
          @fact (g - g0) |> sum |> abs --> 1
        end

        if typeof(l) != NKLandscape
          context("Neutral neighbors should have the same fitness") do
            function test_neutral_neighbors(genotype, landscape)
              nbrs = neutral_neighbors(genotype, landscape)
              score = fitness(genotype, landscape)
              @fact size(nbrs)[1] --> landscape.n
              @fact size(nbrs)[2] --> greater_than_or_equal(0)
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
            @fact size(nbrs)[2] --> greater_than_or_equal(0)
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
          if length(nbrs) > 0
            nbr1 = fittest_neighbor(g, l)
            @fact nbr1 --> nbrs[:,end]
          end
        end

        context("Adaptive walks") do
          function test_adaptive_walk(walk_function, genotype, landscape)
            walk = walk_function(genotype, landscape)
            @fact walk.length --> greater_than_or_equal(0)
            for i = 2:walk.length
              @fact fitness(walk.history_list[:,i], landscape) --> greater_than(fitness(genotype, landscape))
            end
          end

          context("Random adaptive walk should terminate and move uphill") do
            test_adaptive_walk(random_adaptive_walk, g, l)
          end

          context("Greedy adaptive walk should terminate and move uphill") do
            test_adaptive_walk(greedy_adaptive_walk, g, l)
          end

          context("Reluctant adaptive walk should terminate and move uphill") do
            test_adaptive_walk(reluctant_adaptive_walk, g, l)
          end
        end

        context("Neutral walks") do
          function test_neutral_walk(walk_function, genotype, landscape)
            walk = walk_function(genotype, landscape)
            @fact walk.length --> greater_than_or_equal(0)
            f = fitness(genotype, landscape)
            for i = 1:walk.length
              f0 = fitness(walk.history_list[:,i], landscape)
              @fact f0 --> roughly(f)
            end
          end

          context("Neutral walk should terminate and move sideways") do
            test_neutral_walk(neutral_walk, g, l)
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

          for n = [1, 10]
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

          for n = [1, 10]
            rp = rand(Population, l, n)
            test_population_fitnesses(rp, l)
            zp = zeros(Population, l, n)
            test_population_fitnesses(zp, l)
          end
        end
      end

      context("Meta populations") do
        context("Should be the correct size") do
          function test_meta_population_size(p::MetaPopulation, l::Landscape, n, k)
            @fact popsize(p) --> n
            @fact popct(p) --> k
            @fact size(p)[1] --> l.n
            @fact size(p)[2] --> n
            @fact size(p)[3] --> k
          end

          for n = [1, 10]
            for k = [1, 10]
              rp = rand(MetaPopulation, l, n, k)
              test_meta_population_size(rp, l, n, k)
              zp = zeros(MetaPopulation, l, n, k)
              test_meta_population_size(zp, l, n, k)
            end
          end
        end

        context("Should compute fitness") do
          function test_meta_population_fitnesses(p::MetaPopulation, l::Landscape)
            fs = popfits(p, l)
            for ip = 1:popct(p)
              for ig = 1:popsize(p)
                @fact fs[ig,ip] --> fitness(p[:,ig,ip], l)
              end
            end
          end

          for n = [1, 10]
            for k = [1, 10]
              rp = rand(MetaPopulation, l, n, k)
              test_meta_population_fitnesses(rp, l)
              zp = zeros(MetaPopulation, l, n, k)
              test_meta_population_fitnesses(zp, l)
            end
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

        context("Migration") do
          context("Should move genotypes between populations") do
            rp = rand(MetaPopulation, l, n, 4)
            np = migrate(rp, 1.0, n / 2.0 |> ceil |> Int)
            @fact np --> not(rp)
          end
        end

        context("Linear migration") do
          context("Should move genotypes between populations") do
            rp = rand(MetaPopulation, l, n, 4)
            np = linmigrate(rp, 1.0, n / 2.0 |> ceil |> Int)
            @fact np --> not(rp)
          end
        end
      end
    end
  end

  context("Enumeration") do
    l = NKLandscape(2, 0)

    function testlevel(f, b)
      if f < 0.5
        @fact b --> 0
      else
        @fact b --> 1
      end
    end

    function testfitness(f, g)
      fg = fitness(g, l)
      @fact f --> roughly(fg)
    end

    context("lsfits") do
      fits = lsfits(l)

      testfitness(fits[1], [1, 1])
      testfitness(fits[2], [1, 2])
      testfitness(fits[3], [2, 1])
      testfitness(fits[4], [2, 2])
    end

    context("fitlevs") do
      levs = fitlevs(l, 2)

      f1 = fitness([1, 1], l)
      testlevel(f1, levs[1])

      f2 = fitness([1, 2], l)
      testlevel(f2, levs[2])

      f3 = fitness([2, 1], l)
      testlevel(f3, levs[3])

      f4 = fitness([2, 2], l)
      testlevel(f4, levs[4])
    end

    context("levcounts") do
      cts = levcounts(l, 2)

      gfits = [
        fitness([1, 1], l),
        fitness([1, 2], l),
        fitness([2, 1], l),
        fitness([2, 2], l)
      ]
      low = 0
      high = 0
      for gf = gfits
        if gf < 0.5
          low += 1
        else
          high += 1
        end
      end

      @fact cts[1] --> low
      @fact cts[2] --> high
    end

    context("neighbors") do
      la = NKLandscape(2, 0, a=3)
      @fact_throws neighbors(0, la)

      nbrs = neighbors(IntGenotype(0), l)
      @fact length(nbrs) --> 2
      @fact IntGenotype(1) --> anyof(nbrs...)
      @fact IntGenotype(2) --> anyof(nbrs...)
    end
  end
end

include("fastfunc.jl")
include("fastnets.jl")

