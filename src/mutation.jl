import Base.Random: rand

export bwmutate, bwmutate!, bsmutate, bsmutate!

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

@doc """bsmutate(p::Population, mutprob::Float64)

Return a new population with a single locus of each individual in the original
population mutated (bit flipped if `ls.a == 2`) with probability `mutprob`.
"""
function bsmutate(p::Population, mutprob::Float64)
  np = Population(p)
  bsmutate!(np, mutprob)
  return np
end

@doc """bsmutate(mp::MetaPopulation, mutprob::Float64)
"""
function bsmutate(mp::MetaPopulation, mutprob::Float64)
  np = MetaPopulation(mp)
  bsmutate!(np, mutprob)
  return np
end

@doc """bsmutate!(g::Genotype)

Mutate a single locus within the given genotype.
"""
function bsmutate!(g::Genotype)
  mask = AlleleMask(1) << rand(0:(g.landscape.n - 1))
  g.alleles = g.alleles $ mask
end

@doc """bsmutate(p::Population, ls::Landscape)

Mutate a population in-place by mutating up to one locus of each bitstring in
the population.
"""
function bsmutate!(p::Population, mutprob::Float64)
  for g = p.genotypes
    if rand() < mutprob
      bsmutate!(g)
    end
  end
end

@doc """bsmutate!(p::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bsmutate!(mp::MetaPopulation, mutprob::Float64)
  for p = mp.populations
    bsmutate!(p, mutprob)
  end
end

