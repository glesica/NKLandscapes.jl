using Distributions

export neighbors, number_neighbors, random_neighbor, neutral_neighbors, fitter_neighbors, fittest_neighbors, fittest_neighbor

function neighbors(g::Genome, ls::Landscape, muts::Int64)
  # TODO: Implement multiple mutations.
  # TODO: Figure out if `Task`s will be GC'd on a `break`.
  # TODO: Compute neighbors in random order as an option.
  # TODO: Implement all four of Nowak and Krug's methods.
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  if length(g) < muts
    error("gene count is less than mutation count")
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

neighbors(g::Genome, ls::Landscape) = neighbors(g, ls, 1)

function number_neighbors(g::Genome, ls::Landscape, muts::Int64)
  # TODO: Implement multiple mutations
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  if length(g) < muts
    error("gene count is less than mutation count")
  end
  (ls.a - 1) * ls.n
end

number_neighbors(g::Genome, ls::Landscape) = number_neighbors(g, ls, 1)

function random_neighbor(g::Genome, ls::Landscape, muts::Int64)
  # TODO: Implement multi-mutations.
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  if length(g) < muts
    error("gene count is less than mutation count")
  end

  loci = rand(1:length(g))
  allele = sample(setdiff(1:ls.a, [g[loci]]))
  g1 = g[:]
  g1[loci] = allele
  return g1
end

random_neighbor(g::Genome, ls::Landscape) = random_neighbor(g, ls, 1)

# Returns a matrix of all neutral (equal fitness) neighbors as the columns of a
# matrix.
function neutral_neighbors(g::Genome, ls::Landscape, muts::Int64)
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls, muts)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls, muts)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  neutrals = fits .== f0
  nbrs[:,neutrals]
end

neutral_neighbors(g::Genome, ls::Landscape) = neutral_neighbors(g, ls, 1)

# TODO: Test to make sure this works when there are no fitter neighbors.
# Returns a matrix of all fitter neighbors as the columns of a matrix, sorted
# from lowest fitness (left) to highest fitness (right).
function fitter_neighbors(g::Genome, ls::Landscape, muts::Int64)
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls, muts)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  fits = zeros(Float64, count)
  i = 1
  for nbr = neighbors(g, ls, muts)
    nbrs[:,i] = nbr
    fits[i] = fitness(nbr, ls)
    i += 1
  end
  betters = fits .> f0
  nbrs = nbrs[:,betters]
  fits = fits[betters]
  nbrs[:,sortperm(fits)]
end

fitter_neighbors(g::Genome, ls::Landscape) = fitter_neighbors(g, ls, 1)

# Returns the `count` fittest neighbors, even if they are less fit than `g`, as
# the columns of a matrix, sorted from lowest fitness (left) to highest fitness
# (right).
function fittest_neighbors(g::Genome, ls::Landscape, count::Int64, muts::Int64)
  if (count > number_neighbors(g, ls, muts))
    error("count is too large")
  end
  best = zeros(Int64, ls.n, count)
  fits = zeros(Float64, count)
  for nbr = neighbors(g, ls, muts)
    f0, i = findmin(fits)
    f1 = fitness(nbr, ls)
    if f1 > f0
      best[:,i] = nbr
      fits[i] = f1
    end
  end
  best[:,sortperm(fits)]
end

fittest_neighbors(g::Genome, ls::Landscape, count::Int64) = fittest_neighbors(g, ls, count, 1)

function fittest_neighbor(g::Genome, ls::Landscape, muts::Int64)
  nbrs = fittest_neighbors(g, ls, 1, muts)
  return length(nbrs) == 0 ? nbrs : nbrs[:,end]
end

fittest_neighbor(g::Genome, ls::Landscape) = fittest_neighbor(g, ls, 1)

