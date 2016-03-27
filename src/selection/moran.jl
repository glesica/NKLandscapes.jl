import Base.Random: rand

export moransel, moransel!

@doc """moransel(p::Population, iters::Int64)

Conduct Moran selection `iters` number of times. See below for more information
about a Moran process. Basically, this means to first select an individual
within the population to "die", then proportionally selection an individual
to replace it.

https://en.wikipedia.org/wiki/Moran_process
"""
function moransel(p::Population, iters::Int64)
  np = Population(p)
  moransel!(np, iters)
  return np
end

function moransel(mp::MetaPopulation, iters::Int64)
  np = MetaPopulation(mp)
  moransel!(np, iters)
  return np
end

@doc """moransel!(p::Population, iters::Int64)

Conduct Moran selection in-place.
"""
function moransel!(p::Population, iters::Int64)
  fs = popfits(p)
  fmax = maximum(fs)
  n = popsize(p)

  for _ = 1:iters
    indkill = rand(1:n)
    indbirth = 0
    while indbirth == 0
      i = rand(1:n)
      w = fs[i] / fmax
      if rand() < w
        indbirth = i
      end
    end

    p.genotypes[indkill] = Genotype(p.genotypes[indbirth])
  end
end

function moransel!(mp::MetaPopulation, iters::Int64)
  for p = mp.populations
    moransel!(p, iters)
  end
end

