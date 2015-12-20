using Distributions

export all_neighbors, number_neighbors, random_neighbor, neutral_neighbors, fitter_neighbors, fittest_neighbors, fittest_neighbor

function neighbors(g::Genome, ls::Landscape)
  # TODO: Figure out if `Task`s will be GC'd on a `break`.
  # TODO: Implement all four of Nowak and Krug's methods.
  if length(g) != ls.n
    error("genome is of wrong size")
  end
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

function number_neighbors(g::Genome, ls::Landscape)
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  (ls.a - 1) * ls.n
end

function all_neighbors(g::Genome, ls::Landscape)
  count = number_neighbors(g, ls)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  i = 1
  for nbr = neighbors(g, ls)
    nbrs[:,i] = nbr
    i += 1
  end
  return nbrs
end

function random_neighbor(g::Genome, ls::Landscape)
  if length(g) != ls.n
    error("genome is of wrong size")
  end

  loci = rand(1:length(g))
  allele = sample(setdiff(1:ls.a, [g[loci]]))
  g1 = g[:]
  g1[loci] = allele
  return g1
end

# Returns a matrix of all neutral (equal fitness) neighbors as the columns of a
# matrix.
function neutral_neighbors(g::Genome, ls::Landscape)
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
  neutrals = fits .== f0
  return nbrs[:,neutrals]
end

# TODO: Test to make sure this works when there are no fitter neighbors.
# Returns a matrix of all fitter neighbors as the columns of a matrix, sorted
# from lowest fitness (left) to highest fitness (right).
function fitter_neighbors(g::Genome, ls::Landscape)
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
  return nbrs[:,sortperm(fits)]
end

# Returns the `count` fittest neighbors, even if they are less fit than `g`, as
# the columns of a matrix, sorted from lowest fitness (left) to highest fitness
# (right).
function fittest_neighbors(g::Genome, ls::Landscape, count::Int64)
  if (count > number_neighbors(g, ls))
    error("count is too large")
  end
  nbrs = zeros(Int64, ls.n, count)
  fits = zeros(Float64, count)
  for nbr = neighbors(g, ls)
    f0, i = findmin(fits)
    f1 = fitness(nbr, ls)
    if f1 > f0
      nbrs[:,i] = nbr
      fits[i] = f1
    end
  end
  return nbrs[:,sortperm(fits)]
end

function fittest_neighbor(g::Genome, ls::Landscape)
  nbrs = fittest_neighbors(g, ls, 1)
  return nbrs[:,end]
end

