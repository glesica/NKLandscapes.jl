export natural_adaptive_walk, random_adaptive_walk, greedy_adaptive_walk,
  reluctant_adaptive_walk, neutral_walk, neutral_or_fitter_walk,
  fitter_then_neutral_walk, fitness_range_walk

# TODO: Add a tabu walk that looks at its own history, possibly only as far back as when it got to a certain level of fitness, and doesn't revisit nodes.
# TODO: Add an optional limit on walk length.
# TODO: Another tabu strategy is to avoid mutating a given locus more than once every k iterations.
# TODO: Consider optionally doing a depth-first search with backtracking to ensure total exploration of a plateau before declaring it a local optimum.

@doc """A record of the path of an adaptive or random walk.

`strategy` - Can be one of :natural, :random, :greedy, :reluctant,
:fitness_range, :fitter_then_neutral for now.

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
  push!(w.history_set, g)
  w.fitnesses = [w.fitnesses; f]
end

function natural_adaptive_walk(g::Genotype)
end

@doc """An adaptive walk in which a random, fitter neighbor is chosen at each
step.
"""
function random_adaptive_walk(g::Genotype)
  g0 = g
  f0 = fitness(g0)
  w = Walk(:random, g0, f0)
  while true
    nbrs = fitter_neighbors(g0)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,rand(1:end)]
    f0 = fitness(g0)

    add_step!(w, g0, f0)
  end
  return w
end

@doc """An adaptive walk in which the fittest neighbor is chosen at each step.
"""
function greedy_adaptive_walk(g::Genotype)
  g0 = g
  f0 = fitness(g0)
  w = Walk(:greedy, g0, f0)
  while true
    nbrs = fitter_neighbors(g0)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,end]
    f0 = fitness(g0)

    add_step!(w, g0, f0)
  end
  return w
end

@doc """An adaptive walk in which the least fit, fitter neighbor is chosen at
each step.
"""
function reluctant_adaptive_walk(g::Genotype)
  g0 = g
  f0 = fitness(g0)
  w = Walk(:random, g0, f0)
  while true
    nbrs = fitter_neighbors(g0)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,1]
    f0 = fitness(g0)

    add_step!(w, g0, f0)
  end
  return w
end

function generic_neutral_walk(g::Genotype, walktype::Symbol, nbrfunc::Function, maxsteps::Int64)
  g0 = g
  f0 = fitness(g0)
  w = Walk(walktype, g0, f0)

  step = 1
  while true
    nbrs = nbrfunc(g0)
    if length(nbrs) == 0
      break
    end

    # Choose a new genotype from neutral neighbors
    # or terminate if there isn't one we haven't yet
    # visited.
    i = findfirst(j -> !(nbrs[:,j] in w.history_set), 1:size(nbrs)[2])
    if i == 0
      break
    end
    g0 = nbrs[:,i]
    f0 = fitness(g0)
    add_step!(w, g0, f0)

    step += 1
    if maxsteps > 0 && step > maxsteps
      break
    end
  end
  return w
end

@doc """
A purely neutral walk in which a random, equally fit neighbor is chosen at each
step.

Note that this walk is unlikely to generate a path of length greater than one
with the `NKLandscape` because it requires exact fitness equality.

To avoid cycles, this walk will not revisit genotypes. It uses the
`history_set` attribute of the `Walk` type to enforce this condition.

The resulting walk will have no more than `maxsteps` steps.
"""
function neutral_walk(g::Genotype, maxsteps::Int64=0)
  generic_neutral_walk(g, :neutral, neutral_neighbors, maxsteps)
end

@doc """
An approximatly neutral walk in which a random neighbor with fitness in the
closed interval [lb, ub] is chosen at each step.  This walk is suitable for
generating an approximately neutral walk with an `NKLandscape` because it
does not require exact fitness equality.

To avoid cycles, this walk will not revisit genotypes. It uses the
`history_set` attribute of the `Walk` type to enforce this condition.

The resulting walk will have no more than `maxsteps` steps.
"""
function fitness_range_walk(g::Genotype, lb::Float64, ub::Float64, maxsteps::Int64=0)
  nbrfunc = (g) -> fitness_range_neighbors(g, lb, ub)
  generic_neutral_walk(g, :fitness_range, nbrfunc, maxsteps)
end

@doc """
A walk in which a random neighbor is chosen from the set of fitter or equal
neighbors at each step. Note that this walk is equivalent to a random adaptive
walk for an `NKLandscape` because the likelihood of equal fitness is
vanishingly small for an `NKLandsape`.

To avoid cycles, this walk will not revisit genotypes. It uses the
`history_set` attribute of the `Walk` type to enforce this condition.

The resulting walk will have no more than `maxsteps` steps.
"""
function neutral_or_fitter_walk(g::Genotype, maxsteps::Int64=0)
  # TODO: Considering optimizing so history_set is cleared when fitness changes.
  nbrfunc = (g) -> fitter_or_equal_neighbors(g, sort=false)
  generic_neutral_walk(g, :neutral, nbrfunc, maxsteps)
end

@doc """
A walk in which a random fitter neighbor is chosen if there is one. If there is
no fitter neighbor, then a random neutral neighbor is chosen. Note that this
walk is equivalent to a random adaptive walk for an `NKLandscape` because the
likelihood of equal fitness is vanishingly small for an NKLandsape.

To avoid cycles, this walk will not revisit genotypes at the same fitness
level. It uses the `history_set` attribute of the `Walk` type to enforce this
condition.

The resulting walk will have no more than `maxsteps` steps.
"""
function fitter_then_neutral_walk(g::Genotype, maxsteps::Int64=0)
  # TODO: Considering optimizing so history_set is cleared when fitness changes.
  function nbrfunc(g)
    nbrs = fitter_neighbors(g, sort=false)
    if length(nbrs) > 0
      return nbrs
    end
    return neutral_neighbors(g)
  end
  generic_neutral_walk(g, :fitter_then_neutral, nbrfunc, maxsteps)
end

