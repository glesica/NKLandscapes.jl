using FactCheck

import NKLandscapes
const NK = NKLandscapes

srand(1)

facts("migration.jl") do
  # We use complex landscapes here to ensure there are no accidental
  # duplicates.
  landscapes = [
    NK.NKLandscape(96, 1),
    NK.NKqLandscape(96, 1, 2),
    NK.NKpLandscape(96, 1, 0.90),
  ]

  function testmigrate(func::Function, ls::NK.Landscape)
    p = rand(NK.MetaPopulation, ls, 5, 10)
    np = func(p, 1.0, 5)
    for popind = 1:NK.popct(p)
      origpop = p.populations[popind]
      newpop = np.populations[popind]
      sharect = 0
      for genoind = 1:NK.popsize(p)
        origgeno = origpop.genotypes[genoind]
        newgeno = newpop.genotypes[genoind]
        if origgeno == newgeno
          sharect += 1
        end
      end
      @fact sharect --> less_than(NK.popsize(p)) "expected migrations"
    end
  end

  function testmigprob(func::Function, ls::NK.Landscape)
    numpops = 10000
    p = rand(NK.MetaPopulation, ls, 1, numpops)
    np = func(p, 0.25, 1)
    sharect = 0
    for popind = 1:NK.popct(p)
      origgeno = p.populations[popind].genotypes[1]
      newgeno = np.populations[popind].genotypes[1]
      if origgeno == newgeno
        sharect += 1
      end
    end
    sharepct = sharect / numpops
    @fact sharepct --> roughly(0.75, 0.15)
  end

  function testmigct(func::Function, ls::NK.Landscape)
    p = rand(NK.MetaPopulation, ls, 2, 3)
    np = func(p, 1.0, 1)
    for popind = 1:NK.popct(p)
      origpop = p.populations[popind]
      newpop = np.populations[popind]
      sharect = sum(origpop.genotypes .== newpop.genotypes)
      @fact sharect --> 1
    end
  end

  for ls = landscapes

    context("$(ls)") do
      context("migrate(...)") do
        testmigrate(NK.migrate, ls)

        context("should respect migration probability") do
          testmigprob(NK.migrate, ls)
        end

        context("should respect migration count") do
          testmigct(NK.migrate, ls)
        end
      end

      context("linmigrate(...)") do
        testmigrate(NK.linmigrate, ls)

        context("should respect migration probability") do
          testmigprob(NK.linmigrate, ls)
        end

        context("should respect migration count") do
          testmigct(NK.linmigrate, ls)
        end

        context("should only migrate left and right") do
          p = rand(NK.MetaPopulation, ls, 1, 3)
          np = NK.linmigrate(p, 1.0, 1)
          for popind = 1:NK.popct(p)
            newpop = np.populations[popind]
            leftorig = if popind == 1
              p.populations[NK.popct(p)]
            else
              p.populations[popind - 1]
            end
            rightorig = if popind == NK.popct(p)
              p.populations[1]
            else
              p.populations[popind + 1]
            end
            newgeno = newpop.genotypes[1]
            origgenos = [leftorig.genotypes[1], rightorig.genotypes[1]]
            @fact newgeno --> anyof(origgenos...)
          end
        end
      end
    end
  end
end

