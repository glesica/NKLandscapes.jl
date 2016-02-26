import Base.Random: rand, zeros

export MetaPopulation, popsize, popct, popfits

@doc """
Meta populations can be used to simulation spatially structured populations,
such as those that are geographically isolated from one another.

A meta population can be used in all instances where an ordinary population is
valid.
"""
type MetaPopulation
  populations::Vector{Population}
end

@doc """rand(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)

Create a random meta population based on the given landscape. Each population
will consist of `popsize` individuals and there will be `popct` individual
populations.
"""
rand(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64) =
    MetaPopulation([rand(Population, ls, popsize) for _ = 1:popct])

@doc """zeros(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64)

Create a meta population in which all loci "zeroed". Each population will
consist of `popsize` individuals and there will be `popct` individual
populations.
"""
zeros(::Type{MetaPopulation}, ls::Landscape, popsize::Int64, popct::Int64) =
    MetaPopulation([zeros(Population, ls, popsize) for _ = 1:popct])

@doc """popsize(p::MetaPopulation)

Return the number of individuals in each population within the meta population.
"""
popsize(p::MetaPopulation) = length(p.populations[1].genotypes)

@doc """popct(p::MetaPopulation)

Return the number of populations within the meta population.
"""
popct(p::MetaPopulation) = length(p.populations)

@doc """popfits(p::MetaPopulation, ls::Landscape)

Return a matrix of fitness values where the ith row represents the fitness
values of the ith member of each population. Each column represents a
population.
"""
function popfits(p::MetaPopulation)
  fs = zeros(Float64, popsize(p), popct(p))
  for ip = 1:popct(p)
    for ig = 1:popsize(p)
      fs[ig, ip] = fitness(p.populations[ip].genotypes[ig])
    end
  end
  return fs
end

