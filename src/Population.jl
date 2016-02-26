import Base.Random: rand, zeros

export Population, popsize, popfits

typealias Population Vector{Genotype}

@doc """ function rand(::Type{Population}, ls::Landscape, count::Int64)

Returns a random population of genotypes from landscape ls
"""
function rand(::Type{Population}, ls::Landscape, count::Int64)
  Genotype[rand(Genotype,ls) for _=1:count]
end

@doc """ function zeros(::Type{Population}, ls::Landscape, count::Int64)

Returns a population of zero genotypes from landscape ls
"""
function zeros(::Type{Population}, ls::Landscape, count::Int64)
  Genotype[zeros(Genotype,ls) for _=1:count]
end

popsize(p::Population) = length(p)

@doc """function popfits(p::Population)

Returns a vector of the fitnesses of population p
"""
function popfits(p::Population)
  n = popsize(p)
  Float64[ fitness(p[i]) for i=1:n]
end

