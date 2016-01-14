export natural_adaptive_walk, random_adaptive_walk, greedy_adaptive_walk, reluctant_adaptive_walk, neutral_walk

type Walk
  # Can be one of :natural, :random, :greedy, :reluctant for now.
  strategy::Symbol

  # The number of transitions from the original genotype to the ultimate optimum.
  # This will be one less than the width of the history_list matrix.
  length::Int64

  # Array of genotypes found along the path from the original to the optimum.
  # Each genotype is stored in a column.
  history_list::Matrix{Int64}

  # Dictionary of genotypes found along the path from the original to the optimum.
  # Used in plateau search to avoid repeated states
  history_set::Dict{Array{Int64,1},Int64}

  # Vector of fitnesses of the genotype visited along the path to the optimum.
  # The ith fitness corresponds to the ith column in the history_list matrix.
  fitnesses::Vector{Float64}
end

Walk(s::Symbol, g::Genotype, f::Float64) = Walk(s, 0, reshape(g, (length(g), 1)), Dict{Array{Int64,1},Int64}(g=>1), [f])

function add_step!(w::Walk, g::Genotype, f::Float64)
  w.length += 1
  w.history_list = [w.history_list g]
  w.fitnesses = [w.fitnesses; f]
end

function natural_adaptive_walk(g::Genotype, l::Landscape)
end

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
    while i < size(nbrs)[2] && g0 in keys(w.history_set)
      i += 1
      g0 = nbrs[:,i]
    end
    if i >= size(nbrs)[2]
      break
    end
    f0 = fitness(g0, l)
    w.history_set[g0] = 1
    add_step!(w, g0, f0)
  end
  return w
end
