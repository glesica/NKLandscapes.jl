export Basin, basins, basinlists! 
using DataStructures
using Base.Test

@doc """ type Basin

A basin of attraction of a local optimum under an adaptive walk (currently only a greedy adaptive walk).
"""

type Basin
  gtype::Genotype
  count::Int64
  peak_fitness::Float64
end

@doc """ function basins(ls::NKLandscape; fits::Vector{Float64}=zeros(Float64,0))

Returns an instance of `IntDisjointSets` that contains, as disjoint sets, the
basins of attraction (relative to a greedy adaptive walk) present in the given landscape.
The root of the basin of attraction is not necessarily the corresponding local optimum
(but this property holds after basinlists!() is callled).
The optional argument fits is an array of the fitnesses of the genotypes of the landscape.
Should only be applied to NKLandscapes (rather than NKpLandscapes or NKqLandscapes at this time.
"""
function basins(ls::NKLandscape, fits::Vector{Float64}=zeros(Float64,0))
  basin_sets = IntDisjointSets(ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    if length(fits) > 0 
      fg = fits[i+1]
      @assert fg == fitness(g)
    else
      fg = fitness(g)
    end
    fittest = fittest_neighbor(g)
    if length(fits) > 0
      f_fittest = fits[fittest.alleles+1]
      @assert f_fittest == fitness(fittest)
    else
      f_fittest = fitness(fittest)
    end
    if f_fittest <= fg
      continue
    else
      union!(basin_sets,g.alleles+1,fittest.alleles+1)
    end
  end
  return basin_sets
end

#= This version uses a DisjointSets data structure rather than an IntDisjointSets data structure
For the time being, it is "archived" here as an alternative, perhaps more readable, approach.
function basins(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = DisjointSets{IntGenotype}(collect(0:convert(IntGenotype,(ls.a^ls.n-1))))
  for ig::IntGenotype = 0:(ls.a^ls.n - 1)
    fg = fits[ig+1]
    g = itog(ig,ls)
    fittest = fittest_neighbor(g,ls)
    ifittest = gtoi(fittest,ls)
    if fits[ifittest+1] <= fg
      #println(" no union")
      continue
    else
      union!(basin_sets,ig,ifittest)
      #println("    union")
    end
  end
  return basin_sets
end
=#

@doc """ function local_max(g::Genotype)  

Find the fitness of the corresponding local maximum
"""
function local_max(g::Genotype)  
  return greedy_adaptive_walk(g).history_list[end]
end

@doc """function change_root!(s::IntDisjointSets, new_root::Int64)

Changes the root of the set containing "new_root" to be "new_root". 
Does not change the membership of any set.
d.ranks[new_root] may be too large, but this won't change the correctness of any subsequent operation.
"""
function change_root!(s::IntDisjointSets, new_root::UInt128)
  new_root = convert(Int64,new_root)
  old_root::Int64 = find_root(s,new_root)
  s.parents[old_root] = new_root
  s.parents[new_root] = new_root
end

@doc """ function basinlists!(basins::IntDisjointSets,ls::Landscape)

Returns a list of basins of attractions (each of type Basin) sorted by fitness.
Also modifies the IntDisjointSets data structure  s  so that the root of each set
is the corresponding peak (local optimum under greedy adaptive walk).
"""
function basinlists!(basins::IntDisjointSets,ls::NKLandscape)
  c = counter(Int64)
  for i = 0:(ls.a^ls.n-1)
    push!(c,find_root(basins,i+1)-1)  
  end
  basin_list = Basin[]
  for k in keys(c)
    g = local_max(Genotype(k,ls))
    fit_g = fitness(g)
    change_root!(basins,g.alleles+1)
    push!(basin_list,Basin(g,c[k],fit_g))
  end
  sort!(basin_list, by=x -> x.peak_fitness)  # Sort by fitness
  return basin_list
end

@doc """ function ranks(s::IntDisjointSets)

returns the ranks array of the IntDisjointSets data structure s.  Not exported.
"""
function ranks(s::IntDisjointSets)
  return s.internal.ranks
end

@doc """ function parents(s::IntDisjointSets)

returns the parents array of the IntDisjointSets data structure s.  Not exported.
"""
function parents(s::IntDisjointSets)
  return s.internal.parents
end
