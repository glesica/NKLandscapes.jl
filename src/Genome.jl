export Genome, fitness, random_genome

typealias Genome Vector{Int64}

# TODO: Allow random fitnesses to be drawn from other distributions.

contrib(i::Int64, g::Genome, ls::Landscape, update::Function) = begin
  cont_is = ls.links[:,i]
  cont_vs = g[cont_is]
  return get!(update, ls.contribs, [i; cont_vs])
end
contrib(i::Int64, g::Genome, ls::NKLandscape) = contrib(i, g, ls, () -> rand())
contrib(i::Int64, g::Genome, ls::NKqLandscape) = contrib(i, g, ls, () -> rand(0:(ls.q - 1)))
contrib(i::Int64, g::Genome, ls::NKpLandscape) = begin
  update = () -> begin
    if rand() < ls.p
      return 0.0
    else
      return rand()
    end
  end
  return contrib(i, g, ls, update)
end

function fitness(g::Genome, ls::NKqLandscape)
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  (map(1:length(g)) do i
    contrib(i, g, ls)
  end |> mean) / (ls.q - 1)
end

function fitness(g::Genome, ls::Landscape)
  # TODO: Do we include the fake zeros in the sum or just let the fitness range lower?
  if length(g) != ls.n
    error("genome is of wrong size")
  end
  map(1:length(g)) do i
    contrib(i, g, ls)
  end |> mean
end

random_genome(ls::Landscape) = [rand(1:ls.a) for _ = 1:ls.n]

