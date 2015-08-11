export Genome, fitness, random_genome

typealias Genome Vector{Int64}

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

random_genome(ls::Landscape) = [rand(1:ls.a) for _ = 1:ls.n]

