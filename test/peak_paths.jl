# To run stand-alone:  julia -e "include(\"../src/NKLandscapes.jl\"); include(\"peak_paths.jl\")"
using Base.Collections
using DataStructures
using NKLandscapes
using FactCheck

include("../src/basins.jl")
include("../src/peak_paths.jl")

# Needs considerable "clean up", such as deleting debugging statements, improving documentation, deleting debugging functions, etc.

const fit_diff_weight = 1.0

# Utility function 
cv(x) = convert(IntGenotype,x)

type QueueNode
  current::IntGenotype
  previous::Cons
  cost::Float64
end    
QueueNode(cur::Int64, prev::Int64, cost::Float64) = QueueNode( cv(cur), cv(prev), cost )

@doc """

Path holds a list of genotypes (the actual path) and the cost of the path.
To access the first genotype, use  path[1].
To access the last genotype, use  path[end].
"""
type PathWithCost
  path::Array{IntGenotype,1}
  cost::Float64
end

@doc """function test_paths(n,k,num_landscapes)
Finds least-cost paths between all pairs of peaks of "num_landscapes" landscapes.
least-cost paths are computed in both directions as a check for correctness.

Results are printed for each landscape, and summarized for all landscapes.

TODO:  This function combines testing and reporting.  These functions should be separated.
TODO:  Each landscape could be a separate thread.
"""
function test_paths(n,k)
  @assert k > 0
  ls = NKLandscape(2,0)
  bcc = Basin[]
  fts = Array{Float64,1}
  done = false
  while !done
    ls = NKLandscape(n,k)
    fts = lsfits(ls)
    bcc = basincounts(basins(ls,fts),ls)
    if length(bcc) > 1
      done = true
    end
  end
  spw_dict = Dict{Tuple{IntGenotype,IntGenotype},PathWithCost}()
  i = 1
  for i = 1:length(bcc)
    p1 = bcc[i].gtype
    p_set = setdiff!(Set(map(b->b.gtype,bcc)),p1)
    spw = least_cost_paths(p1,p_set,ls,fts)
    for sp in spw   # for each least-cost path
      opposite_sp = get(spw_dict,(sp.path[end],sp.path[1]),Void)
      if opposite_sp != Void   # The opposite direction path has already been computed
        @fact length(sp.path) --> length(opposite_sp.path) "Expected least-cost paths in opposite directions to have the same length"
        @fact sp.cost --> roughly(opposite_sp.cost) "Expected least-cost paths in opposite directions to have the same cost"
      else
        spw_dict[(sp.path[1],sp.path[end])] = sp
      end
    end
  end
end

n = 6
k = 3   # k must be at least 1
test_paths(n,k)

