using Base.Collections
using Base.Test
using DataStructures
using NKLandscapes

@doc """test_basins(ls::Landscape, fits::Vector{Float64})

Tests the `basins` function by checking that the final genotype of a greedy adpative walk 
started from each possible genotype is in the same disjoint set as the starting genotype
and is the root of this disjoint set.
"""
function test_basins(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = basins(ls,fits)
  basin_lists = basinlists!(basin_sets,ls)
  local_max(x::Int64)= greedy_adaptive_walk(Genotype(x,ls)).history_list[end]
  fit_local_max(x) = fitness(local_max(x))
  for i::Int64 = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    greedy_walk = greedy_adaptive_walk(g)
    gopt = greedy_walk.history_list[end]
    @test find_root(basin_sets,i+1) == find_root(basin_sets,gopt.alleles+1)
    @test find_root(basin_sets,i+1) == convert(Int64,gopt.alleles+1)
  end
  # Chekc that the representative of each basin is a local maximum.
  for b in basin_lists
    fit_b = fitness(b.gtype)
    fn = fittest_neighbor(b.gtype)
    fit_fn = fitness(fn)
    @test fit_fn < fit_b
  end
end

ls = NKLandscape(5,3)
fits = lsfits(ls)
test_basins(ls, fits)
