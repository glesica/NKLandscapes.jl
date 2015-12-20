export natural_walk, random_walk, greedy_walk, reluctant_walk

type Walk
  # Can be one of :natural, :random, :greedy, :reluctant for now.
  strategy::Symbol

  # The number of transitions from the original genome to the ultimate optimum.
  # This will be one less than the width of the history matrix.
  length::Int64

  # Array of genomes found along the path from the original to the optimum.
  # Each genome is stored in a column.
  history::Matrix{Int64}

  # Vector of fitnesses of the genomes visited along the path to the optimum.
  # The ith fitness corresponds to the ith column in the history matrix.
  fitnesses::Vector{Float64}
end

Walk(s::Symbol, g::Genome, f::Float64) = Walk(s, 0, reshape(g, (length(g), 1)), [f])

function add_step!(w::Walk, g::Genome, f::Float64)
  w.length += 1
  w.history = [w.history g]
  w.fitnesses = [w.fitnesses; f]
end

function natural_walk(g::Genome, l::Landscape)
end

function random_walk(g::Genome, l::Landscape)
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

function greedy_walk(g::Genome, l::Landscape)
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

function reluctant_walk(g::Genome, l::Landscape)
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
