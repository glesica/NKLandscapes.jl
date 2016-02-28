import Base.Random: rand

export tournsel, tournsel!

@doc """tournsel(p::Population, k::Int64)

Conduct tournament selection on the population and return a new population. For
now, we assume `p=1`, so the best individual in the tournament is selected as a
member of the new population.  Each tournament involves `k` members of the
population.

Note that sampling is done with replacement, so even if `k` is equal to the
population size, there is no guarantee that every individual will participate
in each tournament.

The current implementation is incredibly inefficient and should not be used.

References:

https://en.wikipedia.org/wiki/Tournament_selection
"""
function tournsel(p::Population, k::Int64)
  np = Population(p)
  tournsel!(np, k)
  return np
end

function tournsel(mp::MetaPopulation, k::Int64)
  np = MetaPopulation(mp)
  tournsel!(np, k)
  return np
end

@doc """tournsel!(p::Population, k::Int64)

Conduct tournament selection in-place.
"""
function tournsel!(p::Population, k::Int64)
  n = popsize(p)
  fs = popfits(p)
  selected = zeros(Int64, n)
  for i = 1:n
    # FIXME: This is slow, probably because of indmax over an unsorted list
    contestants = rand(1:n, k)
    winner = fs[contestants] |> indmax
    selected[i] = contestants[winner]
  end
  p.genotypes[:] = p.genotypes[selected]
end

function tournsel!(mp::MetaPopulation, k::Int64)
  for p = mp.populations
    tournsel!(p, k)
  end
end
