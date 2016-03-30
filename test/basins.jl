using Base.Collections
using DataStructures
using NKLandscapes

@doc """test_basins(ls::Landscape, fits::Vector{Float64})

Tests the `basins` function by checking that the final genotype of a greedy adpative walk 
started from each possible genotype is in the same disjoint set as the starting genotype
and is the root of this disjoint set.
"""
function test_basins(ls::NKLandscape, fits::Vector{Float64}=zeros(Float64,0))
  if length(fits) == 0
    fits = lsfits(ls)
  end
  #basin_sets = basins(ls,fits)
  (basin_lists, basin_sets) = basinlist(ls)
  local_max(x::Int64) = greedy_adaptive_walk(Genotype(x, ls)).history_list[end]
  fit_local_max(x) = fitness(local_max(x))
  for i::Int64 = 0:(ls.a^ls.n - 1)
    g = Genotype(i, ls)
    greedy_walk = greedy_adaptive_walk(g)
    gopt = greedy_walk.history_list[end]
    @fact find_root(basin_sets, i + 1) --> find_root(basin_sets, gopt.alleles + 1)
    @fact find_root(basin_sets, i + 1) --> convert(Int64,gopt.alleles + 1)
  end
  # Check that the representative of each basin is a local maximum.
  for b in basin_lists
    fit_b = fitness(b.gtype)
    fn = fittest_neighbor(b.gtype)
    fit_fn = fitness(fn)
    @fact fit_fn --> less_than(fit_b)
  end
end

srand(1)

context("basins.jl") do
  landscapes = [
    NK.NKLandscape(4, 1)
  ]

  for ls = landscapes
    fits = lsfits(ls)

    context("$(ls)") do
      context("basins(...)") do
        test_basins(ls, fits)
      end
    end
  end
end

