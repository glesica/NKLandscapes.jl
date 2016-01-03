import Base.Random: rand, zeros

export Population, popsize, popfits

typealias Population Matrix{Int64}

function rand(::Type{Population}, ls::Landscape, count::Int64)
  p = zeros(Int64, ls.n, count)
  for i = 1:count
    p[:,i] = rand(Genotype, ls)
  end
  p
end

function zeros(::Type{Population}, ls::Landscape, count::Int64)
  p = zeros(Int64, ls.n, count)
  for i = 1:count
    p[:,i] = zeros(Genotype, ls)
  end
  p
end

popsize(p::Population) = size(p)[2]

# Return a vector of fitness values where the ith element in the vector
# corresponds to the fitness of the ith column of the population matrix.
function popfits(p::Population, ls::Landscape)
  n = popsize(p)
  fs = zeros(Float64, n)
  for i = 1:n
    fs[i] = fitness(p[:,i], ls)
  end
  fs
end

