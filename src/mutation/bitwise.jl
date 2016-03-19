import Base.Random: rand
using Distributions
using StatsBase

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

@doc """bwmutate!(g::Genotype, mutprob::Float64)

Mutate the genotype in-place by mutating each locus with probability `mutprob`.
The number of bits mutated is determined by sampling from a Poisson distribution
with mean g.landscape.n*mutprob.  The cases where the number of bits
mutated is 0, 1, or 2 are handled separately for efficiency reasons.
For mutprob >= 0.25, each bit is mutated indendently (which is very inefficient
for small values of mutprob).

Returns the alleles of the modified genotype.
TODO:  Maybe should not return any value?
"""
function bwmutate!(g::Genotype, mutprob::Float64)
  if mutprob <= 0.25  # This is the common case
    num_bits_mutated = rand(Poisson(g.landscape.n*mutprob))
    #println("mutprob: ",mutprob,"  num_bits_mutated: ",num_bits_mutated)
    #println("before: ",g)
    if num_bits_mutated == 0 # handle this common case here for efficiency
      return 
    end
    one_bit_mask = AlleleMask(1)
    if num_bits_mutated == 1
      mask = one_bit_mask << rand(0:g.landscape.n-1)
    elseif num_bits_mutated == 2  
      ind1 = rand(0:(g.landscape.n-1))
      ind2 = rand(0:(g.landscape.n-2))
      if ind1 <= ind2
        ind2 += 1
      end
      mask = (one_bit_mask << ind1) | (one_bit_mask << ind2)
    elseif num_bits_mutated == 3  
      ind1 = rand(0:(g.landscape.n-1))
      ind2 = rand(0:(g.landscape.n-2))
      ind2 = (ind1 <= ind2) ? ind2+1 : ind2
      ind3 = rand(0:(g.landscape.n-3))
      if ind1 < ind2
        ind3 = (ind1 <= ind3) ? ind3+1 : ind3
        ind3 = (ind2 <= ind3) ? ind3+1 : ind3
      else
        ind3 = (ind2 <= ind3) ? ind3+1 : ind3
        ind3 = (ind1 <= ind3) ? ind3+1 : ind3
      end
      mask = (one_bit_mask << ind1) | (one_bit_mask << ind2) | (one_bit_mask << ind3)
    else  # num_bit_mutated > 3
      num_bits_mutated = min(g.landscape.n, num_bits_mutated) 
      bits_mutated = sample(0:(g.landscape.n-1), num_bits_mutated, replace=false)
      mask = AlleleMask(0)
      for b in bits_mutated
        mask $= one_bit_mask << b
      end
    end
  else  # if mutprob > 0.25 # The above is increasing inaccurate for higher values of mutprob.
    # mutate each bit individually
    mask = AlleleMask(0)
    current = AlleleMask(1)
    for _ = 1:g.landscape.n
      if rand() < mutprob
        mask |= current
      end
      current = current << 1
    end
  end
  g.alleles = g.alleles $ mask
  #println("after:  ",g)
  return 
end

@doc """bwmutate!(p::Population, ls::Landscape, mutprob::Float64)

Mutate the population in-place by mutating each locus of each individual in the
population with probability `mutprob`.
"""
function bwmutate!(p::Population, mutprob::Float64)
  for g = p.genotypes
    #println("bw before: ",g)
    bwmutate!(g, mutprob)
    #println("bw after:  ",g)
  end
  #println("bw pop:  ",p)
end

@doc """bwmutate!(mp::MetaPopulation, ls::Landscape, mutprob::Float64)
"""
function bwmutate!(mp::MetaPopulation, mutprob::Float64)
  for p = mp.populations
    bwmutate!(p, mutprob)
  end
end

