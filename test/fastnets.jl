# Command line:   julia -L "NKLandscapes.jl"   
# To run:  
# julia>  using NKLandscapes
# julia>  include("../test/fastnets.jl")
import NKLandscapes
const NK = NKLandscapes
using FactCheck

n = 10
k = 0
a = 2

# Since k == 0, ls has a single minimum fitness neutral net and a single maximum fitness neutral net
# Works for either an NKq landscape or an NK landscape
#ls = NKqLandscape(n,k,a)
ls = NKLandscape(n,k)
fa = fitness_array(ls)
(min_fit,min_ind) = findmin(fa)  
(max_fit,max_ind) = findmax(fa)  
min_ig = min_ind-1    # an integer genotype of minimum fitness
min_g = int_to_genotype(min_ig,ls)  # genotype of minimum fitness
max_ig = max_ind-1    # an integer genotype of maximum fitness
max_g = int_to_genotype(max_ig,ls)  # genotype of maximum fitness
fl = fitness_levels_array(ls,ls.n,fa)   # use n == ls.n fitness levels

facts("NKLandscapes.jl fast tests") do

  context("NK.neighbors(...)") do
    fe_size = size(fitter_or_equal_neighbors(min_g,ls))[2]
    @fact n --> fe_size "Expected number of fitter or equal neighbors to be N = $n"
    nn_size = size(neutral_neighbors(min_g,ls))[2]
    fn_size = size(fitter_neighbors(min_g,ls))[2]
    @fact n --> nn_size + fn_size "Expected number of neutral nbrs + number of fitter nbrs to be N = $n"
    fit_increment = 1.0/ls.n - eps()
    lb = 0.0
    frn_sum = 0
    for i = 0:ls.n-1
      ub = lb + fit_increment
      frn_size = size(fitness_range_neighbors(min_g,ls,lb,ub))[2]
      frn_sum += frn_size
      lb = ub 
    end
    @fact n --> frn_sum "Expected sum of number of fitness range neighbors to be N = $n"
  end

  context("NK.walks(...)") do
    rand_w = random_adaptive_walk(min_g,ls)
    @fact max_fit --> roughly(rand_w.fitnesses[end]) "Expected final fitness of random adaptive walk to be maximum fitness of landscape which is $max_fit"
    greedy_w = greedy_adaptive_walk(min_g,ls)
    @fact max_fit --> roughly(greedy_w.fitnesses[end]) "Expected final fitness of random adaptive walk to be maximum fitness of landscape which is $max_fit"
    reluct_w = reluctant_adaptive_walk(min_g,ls)
    @fact max_fit --> roughly(reluct_w.fitnesses[end]) "Expected final fitness of random adaptive walk to be maximum fitness of landscape which is $max_fit"
    fit_neutral_w = fitter_then_neutral_walk(min_g,ls)
    @fact max_fit --> roughly(fit_neutral_w.fitnesses[end]) "Expected final fitness of random adaptive walk to be maximum fitness of landscape which is $max_fit"
  end

  context("NK.neutral_nets(...)") do
    lnn = sort(list_neutral_nets(ls,fl),lt=(x,y)->x[3]<y[3])  # neutral nets sorted by fitness
    @fact lnn[end][3] --> greater_than(lnn[end-1][3])  "Expected a single neutral net of maximum fitness"
  end
end

