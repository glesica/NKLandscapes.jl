typealias IntGenotype UInt64

@doc """gtoi(g::Genotype, ls::Landscape)

Convert a `Genotype` to an equivalent integer of base `ls.a`.

This operation is limited by the size of a `UInt64`. For example, if `ls.a` is
two, then the genotype must be of length length than or equal to 64.
"""
function gtoi(g::Genotype, ls::Landscape)
  sum::IntGenotype = g[1] - 1
  for i = 2:length(g)
    sum = ls.a * sum + g[i] - 1
  end
  return sum
end

@doc """itog(n::IntGenotype, ls::Landscape)

Convert an integer of base `ls.a` to an equivalent `Genotype`.
"""
function itog(n::IntGenotype, ls::Landscape)
  [parse(IntGenotype, allele, ls.a) + 1 for allele = base(ls.a, n, ls.a)]
end

@doc """fitness(g::IntGenotype, ls::Landscape)

Returns the fitness of the given integer genotype.
"""
function fitness(g::IntGenotype, ls::Landscape)
  fitness(itog(g, ls), ls)
end

@doc """lsfits(ls::Landscape)

Return a vector containing the fitnesses of all possible genotypes for the
given landscape. The i'th element of the vector will correspond to the
fitness of the genotype that, when converted to an integer using `gtoi`, has
the value `i`.
"""
function lsfits(ls::Landscape)
  fits = zeros(Float64, ls.a^ls.n)
  for i = 0:(l.a^l.n - 1)
    fits[i + 1] = fitness(i, ls)
  end
  return fits
end

@doc """levfits(ls::Landscape, ints::Int64, fits::Vector{Float64}, epsilon::Float64)

Return a vector of fitness levels of each genotype possible on the given
landscape. The fitness level is the index of the bin into which it would fall
assuming the interval [0,1] was divided into `ints` bins.
"""
function levfits(ls::Landscape, ints::Int64, fits::Vector{Float64}, epsilon::Float64)
  levs = zeros(Int64, ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    levs[i + 1] = floor(Int64, (fits[i + 1] + epsilon) * ints)
  end
  return levs
end

@doc """levfits(ls::NKqLandscape, ints::Int64, fits::Vector{Float64})
"""
function levfits(ls::NKqLandscape, ints::Int64, fits::Vector{Float64})
  levfits(ls, ints, fits, eps())
end

@doc """levfits(ls::Landscape, ints::Int64, fits::Vector{Float64})
"""
function levfits(ls::NKLandscape, ints::Int64, fits::Vector{Float64})
  levfits(ls, ints, fits, 0.0)
end

@doc """levfits(ls::Landscape, ints::Int64)
"""
function levfits(ls::Landscape, ints::Int64)
  levfits(ls, ints, lsfits(ls))
end

@doc """levcounts(ls::Landscape, ints::Int64, levs::Vector{Int64})

Return a vector such that the i'th element of the vector contains the number of
genotypes that fall into the i'th fitness level bin.
"""
function levcounts(ls::Landscape, ints::Int64, levs::Vector{Int64})
  counts = zeros(Int64, ints)
  for i = 0:(ls.a^ls.n - 1)
    counts[levs[i + 1]] += 1
  end
  return counts
end

@doc """levcounts(ls::Landscape, ints::Int64, fits::Vector{Float64})
"""
function levcounts(ls::Landscape, ints::Int64, fits::Vector{Float64})
  levcounts(ls, ints, levfits(ls, ints, fits))
end

@doc """levcounts(ls::Landscape, ints::Int64)
"""
function levcounts(ls::Landscape, ints::Int64)
  levcounts(ls, ints, lsfits(ls))
end

# assumes landscape l has arity 2, i. e.,  l.a==2.
@doc """neighbors(g::IntGenotype, ls::Landscape)

Return a vector of all one-mutant neighbors as integers.
"""
function neighbors(g::IntGenotype, l::Landscape)
  if ls.a != 2
    error("Only landscapes with a=2 are currently supported by neighbors().")
  end
  return [gtoi(nbr, ls) for nbr = neighbors(vg, ls)]
end

