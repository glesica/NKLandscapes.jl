export neighbors, neutral_neighbors, better_neighbors, best_neighbors, best_neighbor

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

function neutral_neighbors(g::Genome, ls::Landscape, muts::Int64)
  f0 = fitness(g, ls)
  Task(function ()
    for nbr = neighbors(g, ls, muts)
      f1 = fitness(nbr, ls)
      # TODO: Consider doing an approximate comparison here.
      if f1 == f0
        produce(nbr)
      end
    end
  end)
end

neutral_neighbors(g::Genome, ls::Landscape) = neutral_neighbors(g, ls, 1)

function better_neighbors(g::Genome, ls::Landscape, muts::Int64)
  f0 = fitness(g, ls)
  Task(function ()
    for nbr = neighbors(g, ls, muts)
      f1 = fitness(nbr, ls)
      if f1 > f0
        produce(nbr)
      end
    end
  end)
end

better_neighbors(g::Genome, ls::Landscape) = better_neighbors(g, ls, 1)

function best_neighbors(g::Genome, ls::Landscape, count::Int64, muts::Int64)
  best = zeros(Int64, count, ls.n)
  fits = zeros(Float64, count)
  for nbr = neighbors(g, ls, muts)
    f0, i = findmin(fits)
    f1 = fitness(nbr, ls)
    if f1 > f0
      best[i,:] = nbr
      fits[i] = f1
    end
  end
  # Return them using a generator to keep the interface consistent.
  Task(function ()
    for i = sortperm(fits)
      produce(best[i,:] |> vec)
    end
  end)
end

best_neighbors(g::Genome, ls::Landscape, count::Int64) = best_neighbors(g, ls, count, 1)
best_neighbor(g::Genome, ls::Landscape, muts::Int64) = consume(best_neighbors(g, ls, 1, muts))
best_neighbor(g::Genome, ls::Landscape) = consume(best_neighbors(g, ls, 1))

