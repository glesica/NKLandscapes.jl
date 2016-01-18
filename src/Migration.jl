import Base.Random: rand
import Distributions: sample

export migrate, migrate!, linmigrate, linmigrate!

# Migrate individuals between two populations.
function popmigrate!(p::MetaPopulation, srcpopind::Int64,
    destpopind::Int64, migct::Int64)
  psize = popsize(p)
  srcinds = sample(1:psize, migct, replace=false)
  destinds = sample(1:psize, migct, replace=false)
  p[:,destinds,destpopind] = p[:,srcinds,srcpopind]
end

@doc """migrate(p::MetaPopulation, migprob::Float64, migct::Int64)

Migrate individuals between populations in the given meta population. A
population will participate in migration with probability `migprob` and `migct`
individuals will be migrated between populations.

Populations will be paired randomly for migration.

A new meta population will be returned.
"""
function migrate(p::MetaPopulation, migprob::Float64, migct::Int64)
  np = copy(p)
  migrate!(np, migprob, migct)
  return np
end

@doc """migrate!(p::MetaPopulation, migprob::Float64, migct::Int64)

Conduct migration in-place.
"""
function migrate!(p::MetaPopulation, migprob::Float64, migct::Int64)
  pct = popct(p)
  srcinds = sample(1:pct, migprob * pct |> round |> Int, replace=false)
  for srcind = srcinds
    destind = rand(setdiff(1:pct, srcind))
    popmigrate!(p, srcind, destind, migct)
  end
end

@doc """linmigrate(p::MetaPopulation, migprob::Float64, migct::Int64)

Migrate individuals between populations in the given meta population. A
population will participate in migration with probability `migprob` and `migct`
individuals will be migrated between populations.

Populations will only participate in migrations with the immediate neighbors,
at indices `i+1` and `i-1`, with a periodic boundary condition at the ends.

A new meta population will be returned.
"""
function linmigrate(p::MetaPopulation, migprob::Float64, migct::Int64)
  np = copy(p)
  linmigrate!(np, migprob, migct)
  return np
end

@doc """linmigrate!(p::MetaPopulation, migprob::Float64, migct::Int64)

Conduct linear migration in-place.
"""
function linmigrate!(p::MetaPopulation, migprob::Float64, migct::Int64)
  pct = popct(p)
  srcinds = sample(1:pct, migprob * pct |> round |> Int, replace=false)
  for srcind = srcinds
    destind = rand(if srcind == 1
      [pct, srcind + 1]
    elseif srcind == pct
      [srcind - 1, 1]
    else
      [srcind - 1, srcind + 1]
    end)
    popmigrate!(p, srcind, destind, migct)
  end
end

