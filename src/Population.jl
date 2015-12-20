export Population, fitnesses, random_population

typealias Population Matrix{Int64}

function fitnesses(p::Population, ls::Landscape)
  count = size(p)[2]
  fs = zeros(Float64, count)
  for i = 1:count
    fs[i] = fitness(p[:,i], ls)
  end
  fs
end

function random_population(count::Int64, ls::Landscape)
  p = zeros(Int64, ls.n, count)
  for i = 1:count
    p[:,i] = random_genotype(ls)
  end
  p
end

