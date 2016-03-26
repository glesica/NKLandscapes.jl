using FactCheck

import NKLandscapes
const NK = NKLandscapes

const pcount = 4
const psize = 2

srand(1)

facts("migration.jl") do
  landscapes = [
    NK.NKLandscape(4, 1),
    NK.NKqLandscape(4, 1, 2),
    NK.NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscapes

    context("$(ls)") do
      context("migrate(...)") do
        p = rand(MetaPopulation, ls, 5, 10)
        np = migrate(p, 1.0, 1)
        popsdelta = 0
        for i = 1:popct(p)
          for j = 1:popsize(p)
            g0 = p.populations[i].genotypes[j]
            g1 = np.populations[i].genotypes[j]
            if g0 != g1
              popsdelta += 1
            end
          end
        end
        @fact popsdelta --> greater_than_or_equal(1)
      end

      context("linmigrate(...)") do
        p = rand(MetaPopulation, ls, 2, 3)
        np = migrate(p, 1.0, 1)
      end
    end
  end
end

