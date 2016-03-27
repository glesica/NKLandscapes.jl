import Base.Random: rand

export propsel, propsel!

@doc """propsel(p::Population)

Create a new population through proportional selection on the given
population. Employs the algorithm describe in the Lipowski and
Lipowska paper below.

References:

https://en.wikipedia.org/wiki/Fitness_proportionate_selection

A. Lipowski and D. Lipowska, “Roulette-wheel selection via stochastic
acceptance,” Phys. A Stat. Mech. its Appl., vol. 391, no. 6, pp. 2193–2196,
2012.
"""
function propsel(p::Population)
  np = Population(p)
  propsel!(np)
  return np
end

function propsel(mp::MetaPopulation)
  np = MetaPopulation(mp)
  propsel!(np)
  return np
end

@doc """propsel!(p::Population)

Conduct proportional selection in-place.
"""
function propsel!(p::Population)
  fs = popfits(p)
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

  p.genotypes[:] = [Genotype(p.genotypes[selected[i]]) for i = 1:n]
end

function propsel!(mp::MetaPopulation)
  for p = mp.populations
    propsel!(p)
  end
end
