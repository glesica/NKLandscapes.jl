import Base.Random: rand

export tournsel, tournsel!

@doc """tournsel(p::Population, ls::Landscape, k::Int64)

Conduct tournament selection on the population and return a new
population. For now, we assume `p=1`, so the best individual in
the tournament is selected as a member of the new population.
Each tournament involves `k` members of the population.

Note that sampling is done without replacement, so even if
`k` is equal to the population size, there is no guarantee
that every individual will participate in each tournament.

References:

https://en.wikipedia.org/wiki/Tournament_selection
"""
function tournsel(p::Population, ls::Landscape, k::Int64)
  np = copy(p)
  tournsel!(np, ls, k)
  return np
end

@doc """tournsel!(p::Population, ls::Landscape, k::Int64)

Conduct tournament selection in-place.
"""
function tournsel!(p::Population, ls::Landscape, k::Int64)
  n = popsize(p)
  fs = popfits(p, ls)
  selected = zeros(Int64, n)
  for i = 1:n
    # FIXME: This is slow, probably because of indmax over an unsorted list
    contestants = rand(1:n, k)
    winner = fs[contestants] |> indmax
    selected[i] = contestants[winner]
  end
  p[:,:] = p[:,selected]
end
