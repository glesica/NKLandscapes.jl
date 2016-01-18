import Base.Random: rand, zeros

export MetaPopulation, popsize, popct, popfits

@doc """
Meta populations can be used to simulation spatially structured populations,
such as those that are geographically isolated from one another.

A meta population can be used in all instances where an ordinary population is
valid.
"""
typealias MetaPopulation Array{Int64,3}

@doc """rand(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)

Create a random meta population based on the given landscape. Each population
will consist of `popsize` individuals and there will be `popct` individual
populations.
"""
function rand(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)
  p = zeros(Int64, ls.n, popsize, popct)
  for ip = 1:popct
    for ig = 1:popsize
      p[:,ig,ip] = rand(Genotype, ls)
    end
  end
  return p
end

@doc """zeros(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)

Create a meta population in which all loci "zeroed". Each population will
consist of `popsize` individuals and there will be `popct` individual
populations.
"""
function zeros(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)
  p = zeros(Int64, ls.n, popsize, popct)
  for ip = 1:popct
    for ig = 1:popsize
      p[:,ig,ip] = zeros(Genotype, ls)
    end
  end
  return p
end

@doc """popsize(p::MetaPopulation)

Return the number of individuals in each population within the meta population.
"""
popsize(p::MetaPopulation) = size(p)[2]

@doc """popct(p::MetaPopulation)

Return the number of populations within the meta population.
"""
popct(p::MetaPopulation) = size(p)[3]

@doc """popfits(p::MetaPopulation, ls::Landscape)

Return a matrix of fitness values where the ith row represents the fitness
values of the ith member of each population. Each column represents a
population.
"""
function popfits(p::MetaPopulation, ls::Landscape)
  fs = zeros(Float64, popsize(p), popct(p))
  for ip = 1:popct(p)
    for ig = 1:popsize(p)
      fs[ig,ip] = fitness(p[:,ig,ip], ls)
    end
  end
  return fs
end

