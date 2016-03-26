import Base.Random: rand
import Distributions: sample

export migrate, migrate!, linmigrate, linmigrate!

# Migrate individuals between two populations.
function popmigrate!(dstpop::Population, srcpop::Population, migct::Int64)
  dstinds = sample(1:popsize(dstpop), migct, replace=false)
  srcinds = sample(1:popsize(srcpop), migct, replace=false)
  dstpop.genotypes[dstinds] = [Genotype(g) for g = srcpop.genotypes[srcinds]]
end

@doc """migrate(p::MetaPopulation, migprob::Float64, migct::Int64)

Migrate individuals between populations in the given meta population. A
population will participate in an in-migration with probability `migprob` and
`migct` individuals will be migrated between populations.

Populations will be paired randomly for migration.

A new meta population will be returned.
"""
function migrate(p::MetaPopulation, migprob::Float64, migct::Int64)
  np = MetaPopulation(p)
  migrate!(np, migprob, migct)
  return np
end

@doc """migrate!(p::MetaPopulation, migprob::Float64, migct::Int64)

Conduct migration in-place.
"""
# TODO: Make a pristine copy and take migrants from it?
# TODO: Track in-migrants separately and add them at the end?
# TODO: Second option seems more promising...
function migrate!(p::MetaPopulation, migprob::Float64, migct::Int64)
  for dstpop = p.populations
    if rand() >= migprob
      continue
    end
    srcpops = sample(p.populations, 2, replace=false)
    srcpop = if srcpops[1] != dstpop
      srcpops[1]
    else
      srcpops[2]
    end
    popmigrate!(dstpop, srcpop, migct)
  end
end

@doc """linmigrate(p::MetaPopulation, migprob::Float64, migct::Int64)

Migrate individuals between populations in the given meta population. A
population will participate in an in-migration with probability `migprob` and
`migct` individuals will be migrated between populations.

Populations will only participate in migrations with the immediate neighbors,
at indices `i+1` and `i-1`, with a periodic boundary condition at the ends.

A new meta population will be returned.
"""
function linmigrate(p::MetaPopulation, migprob::Float64, migct::Int64)
  np = MetaPopulation(p)
  linmigrate!(np, migprob, migct)
  return np
end

@doc """linmigrate!(p::MetaPopulation, migprob::Float64, migct::Int64)

Conduct linear migration in-place.
"""
function linmigrate!(p::MetaPopulation, migprob::Float64, migct::Int64)
  count = popct(p)
  for i = 1:count
    if rand() >= migprob
      continue
    end
    dstpop = p.populations[i]
    nbrinds = if i == 1
      [count, i + 1]
    elseif i == count
      [i - 1, 1]
    else
      [i - 1, i + 1]
    end
    srcpop = p.populations[rand(nbrinds)]
    popmigrate!(dstpop, srcpop, migct)
  end
end

