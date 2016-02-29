# To run standalone from NKLandscapes.jl directory:  
#  julia -L "src/NKLandscapes.jl" test/fastnets.jl
import NKLandscapes
const NK = NKLandscapes
using FactCheck

n = 3    # n is arbitrary as long as n >= 3, but larger n may give a better test
k = 0    # k must be 0 for these tests to work
q = 2
#a = 3    # a must be 2 for some of the neutral nets functions to work

# Since k == 0, ls has a single minimum fitness neutral net and a single
# maximum fitness neutral net. Works for either an NKq landscape or an NK
# landscape
type LandscapeProperties
  ls::NK.Landscape
  fa::Array{Float64,1}   # Array of fitnesses indexed by integer genotypes
  fl::Array{Int64,1}     # Array of fitness levels indexed by integer genotypes
  min_g::NK.Genotype # A genotype of minimum fitness
end

function LandscapeProperties(ls::NK.Landscape)
  fa = NK.lsfits(ls)
  fl = NK.fitlevs(ls, ls.n, fa)
  min_g = Genotype(indmin(fa) - 1, ls)
  LandscapeProperties(ls, fa, fl, min_g)
end

lsp_nk2 = LandscapeProperties(NK.NKLandscape(n, k, a=2))
lsp_nk3 = LandscapeProperties(NK.NKLandscape(n, k, a=3))
lsp_nkq2 = LandscapeProperties(NK.NKqLandscape(n, k, q, a=2))
lsp_nkq3 = LandscapeProperties(NK.NKqLandscape(n, k, q, a=3))
# TODO:  neutralnets currently fails for a=3.  Debug.
#lsp_list = [lsp_nk2, lsp_nk3, lsp_nkq2, lsp_nkq3]
lsp_list = [lsp_nk2, lsp_nkq2 ]

facts("NKLandscapes.jl fast neighbors, walks, and neutral net tests") do
  context("NK.neighbors(...)") do
    for lsp in lsp_list
      fe_length = length(NK.fitter_neighbors(lsp.min_g,orequal=true))
      #println("fe_length:",fe_length)
      @fact n --> fe_length "Expected number of fitter or equal neighbors to be N = $n"
      nn_length = length(NK.neutral_neighbors(lsp.min_g))
      fn_length = length(NK.fitter_neighbors(lsp.min_g))
      @fact n --> nn_length + fn_length "Expected number of neutral nbrs + number of fitter nbrs to be N = $n"
      fit_increment = 1.0/lsp.ls.n - eps()
      lb = 0.0
      frn_sum = 0
      for i = 0:lsp.ls.n
        ub = lb + fit_increment
        frn_length = length(NK.fitness_range_neighbors(lsp.min_g,lb,ub))
        frn_sum += frn_length
        lb = ub 
      end
      @fact n --> frn_sum "Expected sum of number of fitness range neighbors to be N = $n"
    end
  end

  context("NK.walks(...)") do
    for lsp in lsp_list
      max_fit = maximum(lsp.fa)
      rand_w = NK.random_adaptive_walk(lsp.min_g)
      @fact max_fit --> roughly(fitness(rand_w.history_list[end]))
        "Expected final fitness of random adaptive walk to be maximum fitness of landscape which is $max_fit"
      greedy_w = NK.greedy_adaptive_walk(lsp.min_g)
      @fact max_fit --> roughly(fitness(greedy_w.history_list[end]))
        "Expected final fitness of greedy adaptive walk to be maximum fitness of landscape which is $max_fit"
      reluct_w = NK.reluctant_adaptive_walk(lsp.min_g)
      @fact max_fit --> roughly(fitness(reluct_w.history_list[end]))
        "Expected final fitness of reluctant adaptive walk to be maximum fitness of landscape which is $max_fit"
      fit_neutral_w = NK.fitter_then_neutral_walk(lsp.min_g)
      @fact max_fit --> roughly(fitness(fit_neutral_w.history_list[end]))
        "Expected final fitness of fitter_then_neutral adaptive walk to be maximum fitness of landscape which is $max_fit"
    end
  end

  context("NK.netcounts(...)") do
    cc = 1
    for lsp in lsp_list
      println("cc:",cc)
      cc += 1
      dsets = NK.neutralnets(lsp.ls, lsp.fl)
      println("dsets:",dsets)
      lnn = NK.netcounts(dsets, lsp.fl)
      println("length(lnn):",length(lnn))
      sort!(lnn, by=x -> x[3])
      if length(lnn) > 1
        @fact lnn[end][3] --> greater_than(lnn[end-1][3])  "Expected a single neutral net of maximum fitness"
      else  # All genotypes have the same fitness
        @fact lsp.ls.a^lsp.ls.n --> lnn[end][2]  "Expected a single neutral net with all genotypes"
      end
    end
  end
end

