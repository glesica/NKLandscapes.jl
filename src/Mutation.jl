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

@doc """bsmutate(p::Population, ls::Landscape)

Return a new population with a single locus of each bitstring in the original
population mutated (bit flipped if `ls.a == 2`).
"""
function bsmutate(p::Population, ls::Landscape)
  np = copy(p)
  bsmutate!(np, ls)
  return np
end

@doc """bsmutate(p::Population, ls::Landscape)

Mutate a population in-place by mutating one locus of each bitstring in the
population.
"""
function bsmutate!(p::Population, ls::Landscape)
  for i = 1:popsize(p)
    j = rand(1:ls.n)
    p[j,i] = mutlocus(p[j,i], ls)
  end
end

