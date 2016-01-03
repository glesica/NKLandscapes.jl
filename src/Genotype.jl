import Base.Random: rand, zeros

export Genotype, fitness, random_genotype

typealias Genotype Vector{Int64}

# TODO: Allow random fitnesses to be drawn from other distributions.

contrib(i::Int64, g::Genotype, ls::Landscape, update::Function) = begin
  cont_is = ls.links[:,i]
  cont_vs = g[cont_is]
  return get!(update, ls.contribs, [i; cont_vs])
end
contrib(i::Int64, g::Genotype, ls::NKLandscape) = contrib(i, g, ls, () -> rand())
contrib(i::Int64, g::Genotype, ls::NKqLandscape) = contrib(i, g, ls, () -> rand(0:(ls.q - 1)))
contrib(i::Int64, g::Genotype, ls::NKpLandscape) = begin
  update = () -> begin
    if rand() < ls.p
      return 0.0
    else
      return rand()
    end
  end
  return contrib(i, g, ls, update)
end

function check_genotype_size(g::Genotype, ls::Landscape)
  if length(g) != ls.n
    error("genotype is of wrong size")
  end
end

function fitness(g::Genotype, ls::NKqLandscape)
  check_genotype_size(g, ls)
  (map(1:length(g)) do i
    contrib(i, g, ls)
  end |> mean) / (ls.q - 1)
end

function fitness(g::Genotype, ls::Landscape)
  # TODO: Do we include the fake zeros in the sum or just let the fitness range lower?
  check_genotype_size(g, ls)
  map(1:ls.n) do i
    contrib(i, g, ls)
  end |> mean
end

rand(::Type{Genotype}, ls::Landscape) = [rand(1:ls.a) for _ = 1:ls.n]
zeros(::Type{Genotype}, ls::Landscape) = [1 for _ = 1:ls.n]
