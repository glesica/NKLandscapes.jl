export Genome, fitness, random_genome

typealias Genome Vector{Int64}

function contrib(i::Int64, g::Genome, ls::Landscape)
  cont_is = ls.links[:,i]
  cont_vs = g[cont_is]
  get!(ls.contribs, [i, cont_vs]) do
    if rand() < ls.p
      return 0.0
    else
      # TODO: Allow this to be drawn from an arbitrary distribution.
      return rand()
    end
  end
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

