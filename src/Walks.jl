export natural_adaptive_walk, random_adaptive_walk, greedy_adaptive_walk, reluctant_adaptive_walk, neutral_walk, neutral_or_fitter_walk,
  fitter_then_neutral_walk, fitness_range_walk

@doc """A record of the path of an adaptive or random walk.

`strategy` - Can be one of :natural, :random, :greedy, :reluctant for now.

`length` - The number of transitions from the original genotype to the ultimate
optimum. This will be one less than the width of the `history_list` matrix.

`history_list` - Array of genotypes found along the path from the original to
the optimum. Each genotype is stored in a column.

`history_set` - Dictionary of genotypes found along the path from the original
to the optimum. Used in plateau search to avoid repeated states

`fitnesses` - Vector of fitnesses of the genotype visited along the path to the
optimum. The ith fitness corresponds to the ith column in the `history_list`
matrix.
"""
type Walk
  strategy::Symbol
  length::Int64
  history_list::Matrix{Int64}
  history_set::Set{Genotype}
  fitnesses::Vector{Float64}
end

Walk(s::Symbol, g::Genotype, f::Float64) =
    Walk(s, 0, reshape(g, (length(g), 1)), Set{Genotype}((g,)), [f])

function add_step!(w::Walk, g::Genotype, f::Float64)
  w.length += 1
  w.history_list = [w.history_list g]
  w.fitnesses = [w.fitnesses; f]
end

function natural_adaptive_walk(g::Genotype, l::Landscape)
end

@doc """An adaptive walk in which a random, fitter neighbor is chosen at each
step.
"""
function random_adaptive_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:random, g0, f0)
  while true
    nbrs = fitter_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,rand(1:end)]
    f0 = fitness(g0, l)

    add_step!(w, g0, f0)
  end
  return w
end

@doc """An adaptive walk in which the fittest neighbor is chosen at each step.
"""
function greedy_adaptive_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:greedy, g0, f0)
  while true
    nbrs = fitter_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,end]
    f0 = fitness(g0, l)

    add_step!(w, g0, f0)
  end
  return w
end

@doc """An adaptive walk in which the least fit, fitter neighbor is chosen at
each step.
"""
function reluctant_adaptive_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:random, g0, f0)
  while true
    nbrs = fitter_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,1]
    f0 = fitness(g0, l)

    add_step!(w, g0, f0)
  end
  return w
end

@doc """A purely neutral walk in which a random, equally fit neighbor is chosen
at each step. Note that this walk is unlikely to generate a path of length
greater than one with the `NKLandscape` because it requires exact fitness
equality.

To avoid cycles, this walk will not revisit genotypes. It uses the
`history_set` attribute of the `Walk` type to enforce this condition.
count_limit is an upper bound on the number of steps.
"""
function neutral_walk(g::Genotype, l::Landscape, count_limit=0::Int64)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:neutral, g0, f0)
  count = 0  # step counter
  while true
    nbrs = neutral_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    i = 1
    g0 = nbrs[:,i]
    while i < size(nbrs)[2] && g0 in w.history_set
      i += 1
      g0 = nbrs[:,i]
    end
    if i >= size(nbrs)[2]
      break
    end
    f0 = fitness(g0, l)
    push!(w.history_set, g0)
    add_step!(w, g0, f0)
    count += 1
    if count_limit > 0 && count >= count_limit
      break
    end
  end
  return w
end

@doc """An approximatly neutral walk in which a random neighbor whose fitness
is between lower_bound and upper_bound (inclusive)  is chosen at each step. 
This walk is suitable for generating an approximately neutral walk
with a `NKLandscape` because it requires exact fitness equality.

To avoid cycles, this walk will not revisit genotypes. It uses the
`history_set` attribute of the `Walk` type to enforce this condition.
count_limit is an upper bound on the number of steps.
"""
function fitness_range_walk(g::Genotype, l::Landscape,  lower_bound::Float64, upper_bound::Float64, count_limit::Int64=0)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:neutral, g0, f0)
  count = 0  # step counter
  while true
    nbrs = fitness_range_neighbors(g0, l, lower_bound, upper_bound )
    if length(nbrs) == 0
      break
    end
    i = 1
    g0 = nbrs[:,i]
    while i < size(nbrs)[2] && g0 in w.history_set
      i += 1
      g0 = nbrs[:,i]
    end
    if i >= size(nbrs)[2]
      break
    end
    f0 = fitness(g0, l)
    push!(w.history_set, g0)
    add_step!(w, g0, f0)
    count += 1
    if count_limit > 0 && count >= count_limit
      break
    end
  end
  return w
end

@doc """A walk in which a random neighbor is chosen from the set of
fitter or equal neighbors at each step. Note that this walk is equivalent
to a random adaptive walk for an NKLandscape because the likelihood
of equal fitness is vanishingly small for an NKLandsape.

To avoid cycles, this walk will not revisit genotypes at the same fitness
level. It uses the `history_set` attribute of the `Walk` type to enforce 
this condition.
"""
function neutral_or_fitter_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  #println("f0:",f0,"  g0:",reshape(g0,1,length(g0)))
  w = Walk(:neutral, g0, f0)
  while true
    nbrs = fitter_or_equal_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    i = 1
    g1 = nbrs[:,i]
    f1 = fitness(g1, l)
    while f1 == f0 && i < size(nbrs)[2] && g1 in w.history_set
      i += 1
      g1 = nbrs[:,i]
      f1 = fitness(g1, l)
    end
    #println("f1:",f0,"  g1:",reshape(g1,1,length(g1)))
    if i >= size(nbrs)[2]
      break
    end
    if f1 > f0
      empty!(w.history_set)   
    end
    push!(w.history_set, g1)
    add_step!(w, g1, f1)
    g0 = g1
    f0 = f1
  end
  return w
end

@doc """A walk in which a random fitter neighbor is chosen if there is one.
If there is no fitter neighbor, then a random neutral neighbor is chosen.
Note that this walk is equivalent to a random adaptive walk for an NKLandscape 
because the likelihood of equal fitness is vanishingly small for an NKLandsape.

To avoid cycles, this walk will not revisit genotypes at the same fitness
level. It uses the `history_set` attribute of the `Walk` type to enforce 
this condition.
"""
function fitter_then_neutral_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:neutral, g0, f0)
  while true
    fnbrs = fitter_neighbors(g0, l)
    if length(fnbrs) > 0
      g1 = fnbrs[:,rand(1:end)]
      f1 = fitness(g1, l)
      empty!(w.history_set)
    else
      nnbrs = neutral_neighbors(g0, l)
      if length(nnbrs) == 0
        break
      end
      i = 1
      g1 = nnbrs[:,i]
      f1 = fitness(g1, l)
      while i < size(nnbrs)[2] && g1 in w.history_set
        i += 1
        g1 = nnbrs[:,i]
        f1 = fitness(g1, l)
      end
      if i >= size(nnbrs)[2]
        break
      end
    end
    push!(w.history_set, g1)
    add_step!(w, g1, f1)
    g0 = g1
    f0 = f1
  end
  return w
end

