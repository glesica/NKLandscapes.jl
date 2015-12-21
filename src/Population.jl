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

function popfits(p::Population, ls::Landscape)
  count = popsize(p)
  fs = zeros(Float64, count)
  for i = 1:count
    fs[i] = fitness(p[:,i], ls)
  end
  fs
end

