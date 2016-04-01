#=
We are interested in the difficulty of evolving from one fitness peak to another.  One way that this might happen
is by crossing fitness valleys.  In other words, low fitness genotypes might survive through drift for sufficiently
long to generate mutations in the basin of attraction of a higher fitness peak.  The work of Sergey Gavrilets suggests 
that as the dimension increases, neutral paths between peaks become increasely likely.  However, if these neutral paths
are long, evolution is unlikely to be able to follow them.  So we are interested in paths between peaks which are both
short and approximately neutral.  The non-neutrality of a path might be defined as the sum of the fitness decreases of
steps of the path.  

To be more precise, let g_0, g_1, . . . , g_m  be genotypes of a path.  Let D = f(g_m) - f(g_0) be the net fitness
increase along the path (or decrease if D is negative).  Let F = sum_{i=1}^{m} |f(g_i) - f(g_{i-1})|
= |f(g_1) - f(g_0)| + |f(g_2) - f(g_1)| + . . . + |f(g_m) - f(g_{m-1})|.  It is not hard to see that the sum of
the fitness decreases along the path is  (F - D)/2.

Example:  Let the sequence of fitnesses along the path be 2, 0, 3, 1, 2, 1, 4.  The fitness decreaeses are
(2-0), (3-1), and (2-1), so the non-neutrality is 5.  D = 2 and F = 2 + 3 + 2 + 1 + 1 + 3 = 12.  
(F - D)/2 = (12 - 2)/2 = 10/2 = 5.

The "fit diff" cost F of a path between two genotypes in a landscape is the Hamming distance plus the sum of the
absolute fitness difference between successive genotypes.  

Computes minimum "fit diff" cost paths between all pairs of peaks (local maxima) of NK landscapes.
=#
using Base.Collections
using DataStructures
using NKLandscapes

# Needs considerable "clean up", such as deleting debugging statements, improving documentation, deleting debugging functions, etc.

#const fit_diff_weight = 5.0

type QueueNode
  current::AlleleString
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
  path::Array{AlleleString,1}
  cost::Float64
end

# Utility function for debugging
cv(x) = convert(AlleleString,x)

# For degugging----TODO:  move to test file
function setup(n)
  #include("basins.jl")
  global ls = NKLandscape(n,2)
  global ff = lsfits(ls)
  global bcc = basinlists(basins(ls,ff),ls)
  global p1 = cv(bcc[1].genotype)
  global p2 = cv(bcc[2].genotype)
end

@doc """ function edge_cost(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64})

The cost of a path edge.  
Note the dependence of the cost on the constant "fit_diff_weight".
"""
function edge_cost_exp(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64}, fit_diff_weight::Float64)
  fit_diff = (exp(max(fits[gtype2+1]-fits[gtype1+1],0.0))-exp(0.0))*ls.n*fit_diff_weight
  return 1.0 + fit_diff
end

function edge_cost_sym(gtype1::AlleleString, gtype2::AlleleString, ls::NKLandscape, fits::Vector{Float64}, fit_diff_weight::Float64)
  fit_diff = abs(fits[gtype2+1]-fits[gtype1+1])*ls.n*fit_diff_weight
  return 1.0 + fit_diff
end

@doc """ function least_cost_paths(source::Genotype, destinations::Set, ls::NKLandscape, fits::Vector{Float64}; single_path=false)

Find least cost paths from genotype "source" to the genotypes in the set "destinations" using Dijkstra's algorithm.
"Least cost" is defined by the function "edge_cost_funct".
Returns a list of paths.
If  single_path==true, returns the first path found
 TODO:  Rename variables to be more meaningful.
"""
function least_cost_paths(source::Genotype, destinations::Set{Genotype}, ls::NKLandscape, fits::Vector{Float64}, fit_diff_weight::Float64; 
    edge_cost_funct=edge_cost_exp, single_path=false)
  all_paths = PathWithCost[]
  destination_alleles = Set{AlleleString}([s.alleles for s in destinations])
  closed = Set()
  queue = Base.Collections.PriorityQueue(QueueNode,Float64)
  source_node = QueueNode(source.alleles,cons(source.alleles,nil()),0.0)
  queue[source_node] =0.0
  #println("lcp: queue:",queue)
  while !isempty(queue)
    ng_node = Base.Collections.dequeue!(queue)
    #println("ng_node:",ng_node)
    if ng_node.current in destination_alleles
      setdiff!(destination_alleles, ng_node.current)
      #println("final ng_node:",ng_node.current,"  ",ng_node.cost)
      len = length(ng_node.previous)
      #path = fill(Genotype(0,ls),len)
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
      #println("ng_node:",ng_node.current,"  cost:",ng_node.cost)
      push!(closed, ng_node.current)
      #println("closed:",closed)
      for nbr = neighbors(Genotype(ng_node.current,ls))
        #fit_diff = abs(fitness(nbr,ls)-fitness(ng_node.current,ls))*ls.n*fit_diff_weight
        #heuristic = count_ones( nbr $ peak2 )
        new_cost = edge_cost_funct(nbr.alleles,ng_node.current,ls,fits,fit_diff_weight)
        #println("cur:",ng_node.current,"  nbr:",nbr,"  cur_cost:",ng_node.cost,"   new_cost:",new_cost)
        value = get(queue,nbr,Void)
        if (value == Void) || (value < new_cost)
          nbr_node = QueueNode(nbr.alleles,cons(nbr.alleles,ng_node.previous), new_cost + ng_node.cost)
          #println("nbr:",nbr,"  nbr cost:",nbr_node.cost)
          queue[nbr_node] = new_cost + ng_node.cost
        else
          println("value >:  nbr:",nbr,"  new cost:",nbr_node.cost)
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

function hamming_dist(g1::Genotype,g2::Genotype)
  count_ones(g1.alleles $ g2.alleles)
end

@doc """ function summary_paths(n,k,num_landscapes)

For each of "num_landscapes" landscapes, finds least-cost paths from each peak to all peaks of higher fitness.
Results are printed for each landscape, and summarized for all landscapes.

TODO:  Each landscape could be a separate thread.
"""
function report_paths(n,k,num_landscapes)
  dist_cutoff = 4
  println("fit diff weight: ",fit_diff_weight)
  printlltered_set("dist_cutoff: ",dist_cutoff)
  @printf("   n\t   k\tlen_bcc\tmin_len\tave_len\tave_cst\n")
  sum_ave_length = 0.0
  sum_min_length = 0.0
  sum_ave_cost = 0.0
  sumsq_ave_length = 0.0
  sumsq_min_length = 0.0
  sumsq_ave_cost = 0.0
  for r = 1:num_landscapes
    ls = NKLandscape(n,k)
    fts = lsfits(ls)
    ave_length = 0.0
    min_length = 0.0
    ave_cost = 0.0
    denom = 0
    bcc = sort(basinlists!(basins(ls,fts),ls),lt=(x,y)->x.peak_fitness<y.peak_fitness)
    #println("bcc:",bcc)
    g_set = Set{Genotype}(map(b->b.genotype,bcc))   # List of peak genotypes
    #println("g_set",g_set)
    if length(bcc) > 1
      i = 1
      count = 0
      for i = 1:length(bcc)
        gi = bcc[i].genotype
        #println("gi:",gi)
        setdiff!(g_set,[gi])
        filtered_set = filter(x->hamming_dist(x,gi)<dist_cutoff,g_set)
        #println("i:",i,"  len(g_set):",length(g_set),"  len(f_set):",length(filtered_set))
        spw = least_cost_paths(gi,filtered_set,ls,fts)
        #println("spw:",spw)
        for sp in spw   # for each least-cost path
          ave_length += length(sp.path)-1
          ave_cost += sp.cost
          min_length += count_ones( sp.path[1] $ sp.path[end] )  # Hamming distance between start and end of path
          count += 1
        end
      end
      #denom = length(bcc)*(length(bcc)-1)/2
      denom = count
      ave_length /= denom
      min_length /= denom
      ave_cost /=  denom
      @printf("%4d\t%4d\t%4d\t%6.3f\t%6.3f\t%6.3f\n",n,k,length(bcc),min_length,ave_length,ave_cost)
    else
      # TODO:  Maybe should regenerate the landscape if k is not 0
      println("single peak")
    end
    sum_ave_length += ave_length
    sum_min_length += min_length
    sum_ave_cost += ave_cost
    sumsq_ave_length += ave_length^2
    sumsq_min_length += min_length^2
    sumsq_ave_cost += ave_cost^2
  end
  # TODO: denominator num_landscapes is incorrect if there is a single peak
  std_dev_ave_length = sqrt((sumsq_ave_length - sum_ave_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_min_length = sqrt((sumsq_min_length - sum_min_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_ave_cost = sqrt((sumsq_ave_cost - sum_ave_cost^2/num_landscapes)/(num_landscapes-1))
  @printf("avg_min_length:%6.2f  std_dev_min_length:%6.2f\n",sum_min_length/num_landscapes,std_dev_min_length)
  @printf("avg_ave_length:%6.2f  std_dev_ave_length:%6.2f\n",sum_ave_length/num_landscapes,std_dev_ave_length)
  @printf("avg_ave_cost:%6.2f std_dev_ave_cost:%6.2f\n",sum_ave_cost/num_landscapes,std_dev_ave_cost)
end

function test_paths(ls::NKLandscape, fits::Vector{Float64})
  bcc = sort(basinlists!(basins(ls,fits),ls),lt=(x,y)->x.peak_fitness<y.peak_fitness)
  if length(bcc) > 1
    for i = 1:length(bcc)
      gi = bcc[i].genotype
      gi_set = Set{Genotype}([gi])
      for j = i+1:length(bcc)
        gj = bcc[j].genotype
        gj_set = Set{Genotype}([gj])
        spw_ij = least_cost_paths(gi,gj_set,ls,fits,single_path=true,edge_cost_funct=edge_cost_sym)
        spw_ji = least_cost_paths(gj,gi_set,ls,fits,single_path=true,edge_cost_funct=edge_cost_sym)
        @test  length(spw_ij[1].path) == length(spw_ji[1].path)
        @test_approx_eq spw_ij[1].cost spw_ji[1].cost
      end
    end
  end
end

@doc """print_fitness_diff_matrix(ls::NKLandscape)
  Only for debugging purposes.  Should not be included in any merge with master.
"""
function print_fitness_diff_matrix(ls::NKLandscape)
  print("     ")
  for g = 1:convert(UInt64,(ls.a^ls.n-1))
    @printf("%7d",g)
  end
  println()
  for g = 0:convert(UInt64,(ls.a^ls.n-1))
    @printf("%4d:",g)
    for h = 1:convert(UInt64,(ls.a^ls.n-1))
      fit_diff = abs(fitness(g,ls)-fitness(h,ls))*ls.n*fit_diff_weight
      @printf(" %.4f",fit_diff)
    end
    println()
  end
end
