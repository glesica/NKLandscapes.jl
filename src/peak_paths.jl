#=
We are interested in the difficulty of evolving from one fitness peak to another.  One way that this might happen
is by crossing fitness valleys.  In other words, low fitness genotypes might survive through drift for sufficiently
long to generate mutations in the basin of attraction of a higher fitness peak.  The work of Sergey Gavrilets suggests 
that as the dimension increases, neutral paths between peaks become increasely likely.  However, if these neutral paths
are long, evolution is unlikely to be able to follow them.  So we are interested in paths between peaks which are both
short and approximately neutral.  The non-neutrality of a path might be defined as the sum of the fitness decreases of
steps of the path.  

The cost of a path is the sum of the edge costs of the edges that make up the path.  The edge cost
can be computed in different ways:  two examples are given below.

TODO:  This introductory summary could be improved.
=#
export edge_cost_exp, edge_cost_sym, least_cost_paths, PathWintCost, LandscapePathSummary, summary_landscape_paths,
    report_paths
using Base.Collections
using DataStructures
using NKLandscapes

# Two important parameters:
const fit_diff_weight = 5.0  # Weight on fitness differences as opposed to mutational steps.
const dist_cutoff = 4    # Only paths between pairs of peaks closer to each other than this cutoff are considered

@doc """ type PathWithCost

PathWithCost holds a list of genotypes (the actual path) and the cost of the path.
To access the first genotype, use  path[1].
To access the last genotype, use  path[end].
"""
type PathWithCost
  path::Array{AlleleString,1}
  cost::Float64
end

@doc """ function edge_cost_exp(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64})

The exponential of a cost of a path edge mulitiplied by landscape length and paramter fit_diff_weight
Nonzero only when the fitess differens from gtype1 to gtype2 is negative.
Note the dependence of the cost on the constant "fit_diff_weight".
"""
function edge_cost_exp(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64})
  fit_diff = (exp(max(fits[gtype2+1]-fits[gtype1+1],0.0))-exp(0.0))*ls.n*fit_diff_weight
  return 1.0 + fit_diff
end

@doc """ function edge_cost_sym(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64})

A symmtric (in both directions) edge cost used for testing the correctness of function least_cost_paths().
"""
function edge_cost_sym(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64})
  fit_diff = abs(fits[gtype2+1]-fits[gtype1+1])*ls.n*fit_diff_weight
  return 1.0 + fit_diff
end

@doc """ type QueueNode

Used as a priority queue node in function least_cost_paths()
"""
type QueueNode
  current::AlleleString
  previous::Cons
  cost::Float64
end    
QueueNode(cur::Int64, prev::Int64, cost::Float64) = QueueNode( convert(AlleleString,cur), convert(AlleleString,prev), cost )

@doc """ function least_cost_paths(source::Genotype, destinations::Set, ls::NKLandscape, fits::Vector{Float64}; single_path=false)

Find least cost paths from genotype "source" to the genotypes in the set "destinations" using Dijkstra's algorithm.
"Least cost" is defined by the function "edge_cost_funct".
Returns a list of paths.
If  single_path==true, returns the first path found
 TODO:  Rename variables to be more meaningful.
"""
function least_cost_paths(source::Genotype, destinations::Set{Genotype}, ls::NKLandscape, fits::Vector{Float64}; 
    edge_cost_funct=edge_cost_exp, single_path=false)
  all_paths = PathWithCost[]
  destination_alleles = Set{AlleleString}([s.alleles for s in destinations])
  closed = Set()
  queue = Base.Collections.PriorityQueue(QueueNode,Float64)
  source_node = QueueNode(source.alleles,cons(source.alleles,nil()),0.0)
  queue[source_node] =0.0
  while !isempty(queue)
    ng_node = Base.Collections.dequeue!(queue)
    if ng_node.current in destination_alleles
      setdiff!(destination_alleles, ng_node.current)
      len = length(ng_node.previous)
      path = zeros(AlleleString,len)
      i = len
      for g in ng_node.previous
        #println("g:",g)
        path[i] = g
        i -= 1
      end
      push!(all_paths,PathWithCost(path,ng_node.cost))
      if single_path
        return all_paths
      end
    end
    if !in(ng_node.current, closed)
      push!(closed, ng_node.current)
      for nbr = neighbors(Genotype(ng_node.current,ls))
        #heuristic = count_ones( nbr $ peak2 )  # can be used to convert algorithm to A* for single paths
        new_cost = edge_cost_funct(nbr.alleles,ng_node.current,ls,fits)
        value = get(queue,nbr,Void)
        if (value == Void) || (value < new_cost)
          nbr_node = QueueNode(nbr.alleles,cons(nbr.alleles,ng_node.previous), new_cost + ng_node.cost)
          queue[nbr_node] = new_cost + ng_node.cost
        end
      end
    end
  end
  return all_paths
end

@doc """function path_cost(path::Array{AlleleString,1},ls::NKLandscape)
  An independent computation of the fit-diff cost of a path.
  Mostly for debugging purposes.
"""
function path_cost(path::Array{AlleleString,1},ls::NKLandscape)
  cost = 0.0
  gprev = path[1]
  for gnext in path[2:end]
      hdist = count_ones( gprev $ gnext )
      fit_diff = abs(fitness(gprev,ls)-fitness(gnext,ls))*ls.n*fit_diff_weight
      cost += hdist + fit_diff
      println("hdist:",hdist,"  fit_diff:", ls.n*abs(fitness(gprev,ls)-fitness(gnext,ls)))
      gprev = gnext
  end
  cost
end

@doc """ function hamming_dist(g1::Genotype,g2::Genotype)

Returns Hamming distances between the alleles of g1 and the alleles of g2.
"""
function hamming_dist(g1::Genotype,g2::Genotype)
  count_ones(g1.alleles $ g2.alleles)
end

type LandscapePathSummary
  bcc_length::Int64
  ave_length::Float64
  min_length::Float64
  ave_cost::Float64
  denom::Int64
end

@doc """ function summary_landscape_paths(n::Int64, k::Int64; dist_cutoff::Int64=1000)

Generates an NK landscape with N=n, K=k, and computes all least cost paths for all pairs of
peaks where the first peak is less than the second and where the distance between the peaks
is less than or equal to "dist_cutoff".
Returns a LandscapePathSummary object which summarizes the statistics of these paths.
"""
function summary_landscape_paths(n::Int64, k::Int64; dist_cutoff::Int64=1000)
  if k == 0
    error("k must be greater than or equal to 0 in function summary_landscape_paths")i
  end
  ave_length = 0.0
  min_length = 0.0
  ave_cost = 0.0
  denom = 0
  ls = NKLandscape(n,k)
  fts = lsfits(ls)
  bcc = sort(basinlists!(basins(ls,fts),ls),lt=(x,y)->x.peak_fitness<y.peak_fitness)
  while length(bcc) <= 1
    ls = NKLandscape(n,k)
    fts = lsfits(ls)
    bcc = sort(basinlists!(basins(ls,fts),ls),lt=(x,y)->x.peak_fitness<y.peak_fitness)
  end
  g_set = Set{Genotype}(map(b->b.gtype,bcc))   # List of peak genotypes
  i = 1
  count = 0
  for i = 1:length(bcc)
    gi = bcc[i].gtype
    setdiff!(g_set,[gi])
    filtered_set = filter(x->hamming_dist(x,gi)<dist_cutoff,g_set)
    spw = least_cost_paths(gi,filtered_set,ls,fts)
    for sp in spw   # for each least-cost path
      ave_length += length(sp.path)-1
      ave_cost += sp.cost
      min_length += count_ones( sp.path[1] $ sp.path[end] )  # Hamming distance between start and end of path
      count += 1
    end
  end
  denom = count
  ave_length /= denom
  min_length /= denom
  ave_cost /=  denom
  return LandscapePathSummary(length(bcc), ave_length, min_length, ave_cost, denom)
end

@doc """ function summary_paths(n,k,num_landscapes)

For each of "num_landscapes" landscapes, finds least-cost paths from each peak to all peaks of higher fitness.
Results are printed for each landscape, and summarized for all landscapes.
"""
function report_paths(n,k,num_landscapes)
  println("fit diff weight: ",fit_diff_weight)
  println("dist_cutoff: ",dist_cutoff)
  @printf("   n\t   k\tlen_bcc\tmin_len\tave_len\tave_cst\n")
  sum_ave_length = 0.0
  sum_min_length = 0.0
  sum_ave_cost = 0.0
  sumsq_ave_length = 0.0
  sumsq_min_length = 0.0
  sumsq_ave_cost = 0.0
  path_summary_list = pmap(n_k->summary_landscape_paths(n_k[1], n_k[2], dist_cutoff=5), [(n,k) for _ = 1:num_landscapes])
  for p in path_summary_list
    sum_ave_length += p.ave_length
    sum_min_length += p.min_length
    sum_ave_cost += p.ave_cost
    sumsq_ave_length += p.ave_length^2
    sumsq_min_length += p.min_length^2
    sumsq_ave_cost += p.ave_cost^2
    @printf("%4d\t%4d\t%4d\t%6.3f\t%6.3f\t%6.3f\n",n,k,p.bcc_length,p.min_length,p.ave_length,p.ave_cost)
  end
  std_dev_ave_length = sqrt((sumsq_ave_length - sum_ave_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_min_length = sqrt((sumsq_min_length - sum_min_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_ave_cost = sqrt((sumsq_ave_cost - sum_ave_cost^2/num_landscapes)/(num_landscapes-1))
  @printf("avg_min_length:%6.2f  std_dev_min_length:%6.2f\n",sum_min_length/num_landscapes,std_dev_min_length)
  @printf("avg_ave_length:%6.2f  std_dev_ave_length:%6.2f\n",sum_ave_length/num_landscapes,std_dev_ave_length)
  @printf("avg_ave_cost:%6.2f std_dev_ave_cost:%6.2f\n",sum_ave_cost/num_landscapes,std_dev_ave_cost)
end
