export neighbors, number_neighbors, random_neighbor, neutral_neighbors, better_neighbors, best_neighbors, best_neighbor

function neighbors(g::Genome, ls::Landscape, muts::Int64)
  # TODO: Implement multiple mutations.
  # TODO: Figure out if `Task`s will be GC'd on a `break`.
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
  ls.a * (ls.n - 1)
end

number_neighbors(g::Genome, ls::Landscape) = number_neighbors(g, ls, 1)

function random_neighbor(g::Genome, ls::Landscape, muts::Int64)
  # TODO: Speed this up by jumping directly to the mutant we chose.
  choices = number_neighbors(g, ls, muts)
  k = rand(1:choices)
  i = 1
  for nbr = neighbors(g, ls, muts)
    if i == k
      return nbr
    end
    i += 1
  end
end

random_neighbor(g::Genome, ls::Landscape) = random_neighbor(g, ls, 1)

function neutral_neighbors(g::Genome, ls::Landscape, muts::Int64)
  # TODO: If muts > 1 or N is large this could be a memory burden.
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls, muts)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  neutrals = zeros(Bool, count)
  i = 1
  for nbr = neighbors(g, ls, muts)
    nbrs[:,i] = nbr
    if fitness(nbr, ls) == f0
      neutrals[i] = true
    end
    i += 1
  end
  nbrs[:,neutrals]
end

neutral_neighbors(g::Genome, ls::Landscape) = neutral_neighbors(g, ls, 1)

function better_neighbors(g::Genome, ls::Landscape, muts::Int64)
  # TODO: If muts > 1 or N is large this could be a memory burden.
  f0 = fitness(g, ls)
  count = number_neighbors(g, ls, muts)
  nbrs = zeros(Int64, ls.n, count) # Columns are genotypes
  betters = zeros(Bool, count)
  i = 1
  for nbr = neighbors(g, ls, muts)
    nbrs[:,i] = nbr
    if fitness(nbr, ls) > f0
      betters[i] = true
    end
    i += 1
  end
  nbrs[:,betters]
end

better_neighbors(g::Genome, ls::Landscape) = better_neighbors(g, ls, 1)

function best_neighbors(g::Genome, ls::Landscape, count::Int64, muts::Int64)
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

best_neighbors(g::Genome, ls::Landscape, count::Int64) = best_neighbors(g, ls, count, 1)
best_neighbor(g::Genome, ls::Landscape, muts::Int64) = consume(best_neighbors(g, ls, 1, muts))
best_neighbor(g::Genome, ls::Landscape) = consume(best_neighbors(g, ls, 1))

