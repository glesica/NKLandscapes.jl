import Base.Random: rand

export propsel, propsel!

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
  if fmax == 0
    # all elements have fitness zero
    return
  end

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
