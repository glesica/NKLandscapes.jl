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

function natural_walk(g::Genome, l::Landscape)
end

function random_walk(g::Genome, l::Landscape)
  g0 = g
  f0 = fitness(g0, l)
  w = Walk(:random, 0, reshape(g0, (length(g0), 1)), [f0])
  while true
    nbrs = fitter_neighbors(g0, l)
    if length(nbrs) == 0
      break
    end
    g0 = nbrs[:,rand(1:end)]
    f0 = fitness(g0, l)

    w.length += 1
    w.history = [w.history g0]
    w.fitnesses = [w.fitnesses, f0]
  end
  return w
end

function greedy_walk(g::Genome, l::Landscape)
end

function reluctant_walk(g::Genome, l::Landscape)
end
