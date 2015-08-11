module NK

using Distributions

typealias Genome Vector{Int64}

typealias Population Matrix{Int64}

type Landscape
  n::Int64
  k::Int64
  a::Int64
  links::Matrix{Int64}
  contribs::Array{Float64}
end

function Landscape(n::Int64, k::Int64, p::Float64, a::Int64)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(Int64, n, k + 1)
  for i = 1:n
    links[i,:] = [i, sample(setdiff(1:n, [i]), k, replace=false) |> sort]
  end
  contribs = rand(Float64, n, (ones(Int, k + 1) * a)...)
  # Zero contributions with probability p.
  contribs[rand(Float64, size(contribs)...) .< p] = 0.0
  return Landscape(n, k, a, links, contribs)
end

Landscape(n::Int64, k::Int64) = Landscape(n, k, 0.0, 2)
Landscape(n::Int64, k::Int64, p::Float64) = Landscape(n, k, p, 2)
Landscape(n::Int64, k::Int64, a::Int64) = Landscape(n, k, 0.0, a)

function fitness(g::Genome, ls::Landscape)
  # TODO: Do we include the fake zeros in the sum or just let the fitness range lower?
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  map(1:length(g)) do i
    str_is = ls.links[i,:]
    cont_is = g[str_is]
    ls.contribs[i, cont_is...]
  end |> mean
end

function fitnesses(p::Population, ls::Landscape)
  count = size(p)[2]
  fs = zeros(Float64, count)
  for i = 1:count
    fs[i] = fitness(p[:,i], ls)
  end
  fs
end

random_genome(ls::Landscape) = [rand(1:ls.a) for _ = 1:ls.n]

function random_population(count::Int64, ls::Landscape)
  p = zeros(Int64, ls.n, count)
  for i = 1:count
    p[:,i] = random_genome(ls)
  end
  p
end

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

# TODO: Make a list of mutation algorithms to implement.
end

