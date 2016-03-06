using Base.Collections
using Base.Test
using DataStructures
using NKLandscapes

# Two important parameters:
const fit_diff_weight = 5.0  # Weight on fitness differences as opposed to mutational steps.
const dist_cutoff = 4    # Only paths between pairs of peaks closer to each other than this cutoff are considered

@doc """ function test_paths(ls::NKLandscape, fits::Vector{Float64})

Tests function least_cost_paths() by computing paths between peaks in both
directions using a symmetric edge_cost function, and checking that lengths
and costs are equal.
"""
function test_paths(ls::NKLandscape, fits::Vector{Float64})
  bcc = sort(basinlists!(basins(ls,fits),ls),lt=(x,y)->x.peak_fitness<y.peak_fitness)
  if length(bcc) > 1
    for i = 1:length(bcc)
      gi = bcc[i].gtype
      gi_set = Set{Genotype}([gi])
      for j = i+1:length(bcc)
        gj = bcc[j].gtype
        gj_set = Set{Genotype}([gj])
        spw_ij = least_cost_paths(gi,gj_set,ls,fits,single_path=true,edge_cost_funct=edge_cost_sym)
        spw_ji = least_cost_paths(gj,gi_set,ls,fits,single_path=true,edge_cost_funct=edge_cost_sym)
        @test  length(spw_ij[1].path) == length(spw_ji[1].path)
        @test_approx_eq spw_ij[1].cost spw_ji[1].cost
      end
    end
  end
end

ls = NKLandscape(5,3)
fits = lsfits(ls)
test_paths(ls,fits)
