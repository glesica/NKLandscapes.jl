using DataStructures

export lsfits, fitlevs, levcounts, neighbors, neutralnet, neutralnets, netcounts

@doc """lsfits(ls::Landscape)

Return a vector containing the fitnesses of all possible genotypes for the
given landscape. 
"""
function lsfits(ls::Landscape)
  fits = zeros(Float64, ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    fits[i+1] = fitness(Genotype(i, ls))
  end
  return fits
end

@doc """fitlevs(ls::Landscape, intervals::Int64, fits::Vector{Float64}, epsilon::Float64)

Return a vector of fitness levels of each genotype possible on the given
landscape. The fitness level is the index of the bin into which it would fall
assuming the interval [0,1] was divided into `intervals` bins.
"""
function fitlevs(ls::Landscape, intervals::Int64, fits::Vector{Float64}, epsilon::Float64)
  levels = zeros(Int64, ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    levels[i + 1] = floor(Int64, (fits[i + 1] + epsilon) * intervals)
  end
  return levels
end

@doc """fitlevs(ls::NKqLandscape, intervals::Int64, fits::Vector{Float64})
"""
function fitlevs(ls::NKqLandscape, intervals::Int64, fits::Vector{Float64})
  fitlevs(ls, intervals, fits, eps())
end

@doc """fitlevs(ls::Landscape, intervals::Int64, fits::Vector{Float64})
"""
function fitlevs(ls::NKLandscape, intervals::Int64, fits::Vector{Float64})
  fitlevs(ls, intervals, fits, 0.0)
end

@doc """fitlevs(ls::Landscape, intervals::Int64)
"""
function fitlevs(ls::Landscape, intervals::Int64)
  fitlevs(ls, intervals, lsfits(ls))
end

@doc """levcounts(ls::Landscape, intervals::Int64, levels::Vector{Int64})

Return a vector such that the i'th element of the vector contains the number of
genotypes that fall into the i'th fitness level bin.
"""
function levcounts(ls::Landscape, intervals::Int64, levels::Vector{Int64})
  counts = zeros(Int64, intervals)
  for i = 0:(ls.a^ls.n - 1)
    counts[levels[i + 1] + 1] += 1
  end
  return counts
end

@doc """levcounts(ls::Landscape, intervals::Int64, fits::Vector{Float64})
"""
function levcounts(ls::Landscape, intervals::Int64, fits::Vector{Float64})
  levcounts(ls, intervals, fitlevs(ls, intervals, fits))
end

@doc """levcounts(ls::Landscape, intervals::Int64)
"""
function levcounts(ls::Landscape, intervals::Int64)
  levcounts(ls, intervals, lsfits(ls))
end

@doc """ function neutralnet(g::Genotype, levels::Vector{Int64})

Returns the alleles of the genotypes of the connected neutral network of which the given genotype
is a part as a vector of `AlleleString`s. Two genotypes are considered to be
neighbors if they have the same fitness level. Two genotypes are considered to
be connected if they are neighbors (they differ at one locus).

`levels` is a vector of fitness levels, or bins, as returned by the `fitlevs`
function.
"""
function neutralnet(g::Genotype, levels::Vector{Int64})
  f0 = levels[g.alleles + 1]
  closed = Set{AlleleString}()
  stack = Stack(Genotype)
  push!(stack, g)
  while !isempty(stack)
    ng = pop!(stack)
    if !in(ng.alleles, closed)
      push!(closed, ng.alleles)
      for nbr = filter(nbr -> levels[nbr.alleles + 1] == f0, all_neighbors(ng))
        push!(stack, nbr)
      end
    end
  end
  return collect(closed)
end

@doc """neutralnets(ls::Landscape, levels::Vector{Int64})

Returns an instance of `IntDisjointSet` that contains, as disjoint sets, the
connected neutral networks present in the given landscape.
"""
function neutralnets(ls::Landscape, levels::Vector{Int64})
  nets = IntDisjointSets(ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    if (find_root(nets, i + 1) != i + 1) || (nets.ranks[i + 1] > 0)
      continue   # g has already been accounted for
    end
    net = neutralnet(g, levels)
    if length(net) > 0
      push!(net, i)  # add i to net
    else
      net = AlleleString[i]   # Array containing only i
    end
    for j = 2:length(net)
      union!(nets, net[j - 1] + 1, net[j] + 1)
    end
  end
  return nets
end

@doc """netcounts(nets::IntDisjointSets, levels::Vector{Int64})

Returns a vector of triples, each of which corresponds to a connected neutral
network within the given `IntDisjointSets` instance, which will most likely
come from a call to `neutralnets(...)`. Each triple consists of the following
data:

  1. Integer genotype representative of the network
  2. Number of genotypes in the network
  3. Integer fitness level of the network

The resulting vector is sorted by the number of genotypes in each network.

`levels` is the fitness levels vector that was provided to `neutralnets`.
"""
function netcounts(nets::IntDisjointSets, levels::Vector{Int64})
  c = counter(nets.parents)
  nn_list = map(x -> (x, c[x], levels[x]), keys(c))
  sort!(nn_list, by=x -> x[2]) # Sort by network size
  return nn_list
end

