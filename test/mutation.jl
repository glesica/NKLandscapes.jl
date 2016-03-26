pcount = 4
psize = 2

srand(1)

context("mutation.jl") do
  landscapes = [
    NK.NKLandscape(4, 1),
    NK.NKqLandscape(4, 1, 2),
    NK.NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscapes

    context("$(ls)") do
      context("bwmutate(...)") do
        p = rand(NK.Population, ls, psize)
        np = NK.bwmutate(p, 1.0)
        @fact NK.popsize(np) --> psize
        @fact np --> not(p)
        for i = 1:psize
          @fact np.genotypes[i] --> not(p.genotypes[i])
        end

        mp = rand(NK.MetaPopulation, ls, psize, pcount)
        nmp = NK.bwmutate(mp, 1.0)
        @fact NK.popct(nmp) --> pcount
        @fact NK.popsize(nmp) --> psize
        @fact nmp --> not(mp)
        for i = 1:pcount
          for j = 1:psize
            @fact nmp.populations[i].genotypes[j] --> not(mp.populations[i].genotypes[j])
          end
        end
      end

      context("bsmutate(...)") do
        p = rand(NK.Population, ls, psize)
        np = NK.bsmutate(p, 1.0)
        @fact NK.popsize(np) --> psize
        @fact np --> not(p)
        for i = 1:psize
          @fact np.genotypes[i] --> not(p.genotypes[i])
        end

        mp = rand(NK.MetaPopulation, ls, psize, pcount)
        nmp = NK.bsmutate(mp, 1.0)
        @fact NK.popct(nmp) --> pcount
        @fact NK.popsize(nmp) --> psize
        @fact nmp --> not(mp)
        for i = 1:pcount
          for j = 1:psize
            @fact nmp.populations[i].genotypes[j] --> not(mp.populations[i].genotypes[j])
          end
        end
      end

      context("bsmutate!(::Genotype)") do
        g = rand(NK.Genotype, ls)
        og = deepcopy(g)
        NK.bsmutate!(g)
        @fact g.alleles --> not(og.alleles)
      end
    end
  end
end
