import Base.Random: rand, zeros 

export Population, popsize, popfits, constant

@doc """A group of genotypes.
"""
type Population
  genotypes::Vector{Genotype}
end

Population(p::Population) = Population([Genotype(g) for g = p.genotypes])

@doc """rand(::Type{Population}, ls::Landscape, popsize::Int64)

Returns a random population of genotypes from landscape ls
"""
rand(::Type{Population}, ls::Landscape, popsize::Int64) =
    Population([rand(Genotype, ls) for _ = 1:popsize])

@doc """zeros(::Type{Population}, ls::Landscape, popsize::Int64)

Returns a population of zero genotypes from landscape ls
"""
zeros(::Type{Population}, ls::Landscape, popsize::Int64) =
    Population([zeros(Genotype, ls) for _ = 1:popsize])

@doc """ constant(::Type{Population}, gtype::Genotype, popsize::Int64) =

Returns  population initalized to popsize copies of genotype gtype 
"""
constant(::Type{Population}, gtype::Genotype, popsize::Int64) =
    Population([gtype for _ = 1:popsize])

@doc """popsize(p::Population)

Returns the size of the population.
"""
popsize(p::Population) = length(p.genotypes)

@doc """popfits(p::Population)

Returns a vector of the fitnesses of population p
"""
popfits(p::Population) = [fitness(g) for g = p.genotypes]

