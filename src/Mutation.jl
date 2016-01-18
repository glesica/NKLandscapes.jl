import Base.Random: rand

export bwmutate, bwmutate!, bsmutate, bsmutate!

mutlocus(v::Int64, ls::Landscape) = setdiff(1:ls.a, v) |> rand

@doc """bwmutate(p::Population, ls::Landscape, mutprob::Float64)

Return a new population such that each locus of the original
population has been mutated with probability `mutprob`.
"""
function bwmutate(p::Population, ls::Landscape, mutprob::Float64)
  np = copy(p)
  bwmutate!(np, ls, mutprob)
  return np
end

@doc """bwmutate(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bwmutate(p::MetaPopulation, ls::Landscape, mutprob::Float64)
  np = copy(p)
  bwmutate!(np, ls, mutprob)
  return np
end

@doc """bwmutate!(p::Population, ls::Landscape, mutprob::Float64)

Mutate the population in-place by mutating each locus with probability
`mutprob`.
"""
function bwmutate!(p::Population, ls::Landscape, mutprob::Float64)
  # Use the fact that matrices can be iterated and indexed directly.
  for i = 1:length(p)
    if rand() < mutprob
      p[i] = mutlocus(p[i], ls)
    end
  end
end

@doc """bwmutate!(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bwmutate!(p::MetaPopulation, ls::Landscape, mutprob::Float64)
  for ip = 1:popct(p)
    bwmutate!(p[:,:,ip], ls, mutprob)
  end
end

@doc """bsmutate(p::Population, ls::Landscape)

Return a new population with a single locus of each bitstring in the original
population mutated (bit flipped if `ls.a == 2`) with probability `mutprob`.
"""
function bsmutate(p::Population, ls::Landscape, mutprob::Float64)
  np = copy(p)
  bsmutate!(np, ls, mutprob)
  return np
end

@doc """bsmutate(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bsmutate(p::MetaPopulation, ls::Landscape, mutprob::Float64)
  np = copy(p)
  bsmutate!(np, ls, mutprob)
  return np
end

@doc """bsmutate(p::Population, ls::Landscape)

Mutate a population in-place by mutating up to one locus of each bitstring in
the population.
"""
function bsmutate!(p::Population, ls::Landscape, mutprob::Float64)
  for i = 1:popsize(p)
    if rand() < mutprob
      j = rand(1:ls.n)
      p[j,i] = mutlocus(p[j,i], ls)
    end
  end
end

@doc """bsmutate!(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bsmutate!(p::MetaPopulation, ls::Landscape, mutprob::Float64)
  for ip = 1:popct(p)
    bsmutate!(p[:,:,ip], ls, mutprob)
  end
end

