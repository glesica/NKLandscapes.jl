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

  for ls = landscapes

    context("$(ls)") do
      context("migrate(...)") do
        p = rand(NK.MetaPopulation, ls, 5, 10)
        np = NK.migrate(p, 1.0, 5)
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

        context("should respect migration probability") do
          p = rand(NK.MetaPopulation, ls, 1, 1000)
          np = NK.migrate(p, 0.25, 1)
          sharect = 0
          for popind = 1:NK.popct(p)
            origgeno = p.populations[popind].genotypes[1]
            newgeno = np.populations[popind].genotypes[1]
            if origgeno == newgeno
              sharect += 1
            end
          end
          sharepct = sharect / 1000
          @fact sharepct --> roughly(0.75, 0.25)
        end

        context("should respect migration count") do
          p = rand(NK.MetaPopulation, ls, 2, 2)
          np = NK.migrate(p, 1.0, 1)
          for popind = 1:NK.popct(p)
            origpop = p.populations[popind]
            newpop = np.populations[popind]
            sharect = sum(origpop.genotypes .== newpop.genotypes)
            @fact sharect --> 1
          end
        end
      end

      context("linmigrate(...)") do
        p = rand(NK.MetaPopulation, ls, 5, 10)
        np = NK.migrate(p, 1.0, 1)

        context("should respect migration probability") do
          p = rand(NK.MetaPopulation, ls, 1, 1000)
          np = NK.linmigrate(p, 0.25, 1)
          sharect = 0
          for popind = 1:NK.popct(p)
            origgeno = p.populations[popind].genotypes[1]
            newgeno = np.populations[popind].genotypes[1]
            if origgeno == newgeno
              sharect += 1
            end
          end
          sharepct = sharect / 1000
          @fact sharepct --> roughly(0.75, 0.25)
        end

        context("should respect migration count") do
          p = rand(NK.MetaPopulation, ls, 2, 2)
          np = NK.linmigrate(p, 1.0, 1)
          for popind = 1:NK.popct(p)
            origpop = p.populations[popind]
            newpop = np.populations[popind]
            sharect = sum(origpop.genotypes .== newpop.genotypes)
            @fact sharect --> 1
          end
        end
      end
    end
  end
end

