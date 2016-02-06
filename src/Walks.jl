export natural_adaptive_walk, random_adaptive_walk, greedy_adaptive_walk, reluctant_adaptive_walk, neutral_walk

# TODO: Add a tabu walk that looks at its own history, possibly only as far back as when it got to a certain level of fitness, and doesn't revisit nodes.
# TODO: Add an optional limit on walk length.
# TODO: Another tabu strategy is to avoid mutating a given locus more than once every k iterations.
# TODO: Consider optionally doing a depth-first search with backtracking to ensure total exploration of a plateau before declaring it a local optimum.

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
"""
function neutral_walk(g::Genotype, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:neutral, g0, f0)
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
  end
  return w
end

