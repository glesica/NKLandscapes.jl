import Base.Random: rand

export bwmutate, bwmutate!

@doc """bwmutate(p::Population, ls::Landscape, mutprob::Float64)

Return a new population such that each locus of each individual in the original
population has been mutated with probability `mutprob`.
"""
function bwmutate(p::Population, mutprob::Float64)
  np = Population(p)
  bwmutate!(np, mutprob)
  return np
end

@doc """bwmutate(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bwmutate(mp::MetaPopulation, mutprob::Float64)
  np = MetaPopulation(mp)
  bwmutate!(np, mutprob)
  return np
end

@doc """bwmutate(g::Genotype, mutprob::Float64)
"""
function bwmutate(g::Genotype, mutprob::Float64)
  ng = Genotype(g)
  bwmutate!(ng, mutprob)
  return ng
end

@doc """bwmutate!(g::Genotype, mutprob::Float64)

Mutate the genotype in-place by mutating each locus with probability `mutprob`.
"""
function bwmutate!(g::Genotype, mutprob::Float64)
  mask = AlleleMask(0)
  current = AlleleMask(1)
  for _ = 1:g.landscape.n
    if rand() < mutprob
      mask = mask | current
    end
    current = current << 1
  end
  g.alleles = g.alleles $ mask
end

@doc """bwmutate!(p::Population, ls::Landscape, mutprob::Float64)

Mutate the population in-place by mutating each locus of each individual in the
population with probability `mutprob`.
"""
function bwmutate!(p::Population, mutprob::Float64)
  for g = p.genotypes
    bwmutate!(g, mutprob)
  end
end

@doc """bwmutate!(mp::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bwmutate!(mp::MetaPopulation, mutprob::Float64)
  for p = mp.populations
    bwmutate!(p, mutprob)
  end
end

