import Base.Random: rand, zeros

export Population, popsize, popfits, all_genotypes

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

# Returns a population of all of the possible genotypes for the landscape.
# Requires that the arity ls.a be 2.
# This creates a matrix whose number of columns is exponential in the landscape length ls.n,
#   so it should only run for moderate values of ls.n.
function all_genotypes(::Type{Population}, ls::Landscape )
    @assert  ls.a == 2 "Expected ls.a == 2.  The landscape arity must be 2"
    p = zeros(Int64, ls.n, 2^ls.n)
    for i = 0:(2^ls.n-1)
        p[:,i+1] = int_to_genotype(ls.n,i)
    end
    return p
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


# Convert Int64 to a genotype
# len is the length of the desired genotype
# i is the Int64 to be converted
function int_to_genotype( N, i )
    A = []
    for j = 1:N
        #println("i:",i)
        if i & 1 != 0
            push!(A,2)
        else
            push!(A,1)
        end
        i >>>= 1
    end
    return A
end
