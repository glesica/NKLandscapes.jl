import Base.Random: rand

export bsmutate, bsmutate!

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

