import Base.Random: rand

export propsel, propsel!

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

@doc """propsel(p::Population, ls::Landscape)

Create a new population through proportional selection on the given
population. Employs the algorithm describe in the Lipowski and
Lipowska paper below.

References:

https://en.wikipedia.org/wiki/Fitness_proportionate_selection

A. Lipowski and D. Lipowska, “Roulette-wheel selection via stochastic
acceptance,” Phys. A Stat. Mech. its Appl., vol. 391, no. 6, pp. 2193–2196,
2012.
"""
function propsel(p::Population, ls::Landscape)
  np = copy(p)
  propsel!(np, ls)
  return np
end

@doc """propsel!(p::Population, ls::Landscape)

Conduct proportional selection in-place.
"""
function propsel!(p::Population, ls::Landscape)
  fs = popfits(p, ls)
  fmax = maximum(fs)

  n = popsize(p)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = fs[i] / fmax
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end

  p[:,:] = p[:,selected]
end
