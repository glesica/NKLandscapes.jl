#using Distributions

export all_neighbors, number_neighbors, random_neighbor, neutral_neighbors, fitter_neighbors, fittest_neighbors, fittest_neighbor,
  fitter_or_equal_neighbors, fitness_range_neighbors

function neighbors(g::Genotype, ls::Landscape)
  # TODO: Figure out if `Task`s will be GC'd on a `break`.
  # TODO: Implement all four of Nowak and Krug's methods.
  check_genotype_size(g, ls)
  Task(function ()
    for i = 1:length(g)
      for a = 1:ls.a
        gp = g[:]
        if gp[i] == a
          continue
        end
        gp[i] = a
        produce(gp)
      end
    end
  end)
end

function number_neighbors(g::Genotype, ls::Landscape)
  check_genotype_size(g, ls)
  (ls.a - 1) * ls.n
end

@doc """
Returns all neighbors as the columns of a matrix.
"""
function all_neighbors(g::Genotype, ls::Landscape)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    i += 1
  end
  return nbrs
end

@doc """
Chooses and returns a randomly selected neighbor as a column vector.
"""
function random_neighbor(g::Genotype, ls::Landscape)
  check_genotype_size(g, ls)

  loci = rand(1:length(g))
  allele = sample(setdiff(1:ls.a, [g[loci]]))
  g1 = g[:]
  g1[loci] = allele
  return g1
end

@doc """
Returns all neutral (equal fitness) neighbors as the columns of a matrix. The
neighbors are in random order.
"""
function neutral_neighbors(g::Genotype, ls::Landscape)
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  neutrals = filter(i -> isapprox(fits[i], f0), 1:count) |> shuffle
  return nbrs[:,neutrals]
end

@doc """
Returns all fitter neighbors as the columns of a matrix, sorted
from lowest fitness (left) to highest fitness (right) unless `sort` is set to
`false`, in which case they are randomized.
"""
function fitter_neighbors(g::Genotype, ls::Landscape, sort::Bool=true)
  # TODO: Test to make sure this works when there are no fitter neighbors.
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  betters = fits .> f0
  nbrs = nbrs[:,betters]
  fits = fits[betters]
  return if sort
    nbrs[:,sortperm(fits)]
  else
    nbrs[:,shuffle(fits)]
  end
end

@doc """
Returns all fitter or equal neighbors as the columns of a matrix, sorted from
lowest fitness (left) to highest fitness (right) unless `sort` is set to
`false`, in which case they are randomized.
"""
function fitter_or_equal_neighbors(g::Genotype, ls::Landscape, sort::Bool=true)
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  betters_neutrals = filter(i -> fits[i] >= f0 || isapprox(fits[i], f0), 1:count)
  nbrs = nbrs[:,betters_neutrals]
  fits = fits[betters_neutrals]
  return if sort
    nbrs[:,sortperm(fits)]
  else
    nbrs[:,shuffle(fits)]
  end
end

@doc """
Returns all neighbors with fitness in the closed interval [lb, ub]
as columns of a matrix, sorted from lowest fitness (left) to highest fitness
(right).
"""
function fitness_range_neighbors(g::Genotype, ls::Landscape, lb::Float64, ub::Float64)
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  selected = filter(i -> lb <= fits[i] <= ub, 1:count)
  nbrs = nbrs[:,selected]
  fits = fits[selected]
  return nbrs[:,sortperm(fits)]
end

@doc """
Returns the `count` fittest neighbors, even if they are less fit than `g`, as
the columns of a matrix, sorted from lowest fitness (left) to highest fitness
(right).
"""
function fittest_neighbors(g::Genotype, ls::Landscape, count::Int64)
  if (count > number_neighbors(g, ls))
    error("count is too large")
  end
  nbrs = zeros(Int64, ls.n, count)
  fits = zeros(Float64, count)
  for nbr = neighbors(g, ls)
    f0, i = findmin(fits)
    f1 = fitness(nbr, ls)
    if f1 >= f0
      nbrs[:,i] = nbr
      fits[i] = f1
    end
  end
  return nbrs[:,sortperm(fits)]
end

@doc """
Returns the single fittest neighbor as a column vector.
"""
function fittest_neighbor(g::Genotype, ls::Landscape)
  nbrs = fittest_neighbors(g, ls, 1)
  return nbrs[:,end]
end

