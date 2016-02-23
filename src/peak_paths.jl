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


of the deleterious mutations.
The "fit diff" cost F of a path between two genotypes in a landscape is the Hamming distance plus the sum of the
absolute fitness difference between successive genotypes.  The 
Computes minimum "fit diff" cost paths between all pairs of peaks (local maxima) of NK landscapes.
=#
using Base.Collections
using DataStructures
using NKLandscapes

const fit_diff_weight = 1.0

type QueueNode
  current::IntGenotype
  previous::Cons
  cost::Float64
end    
#QueueNode(cur::UInt64, prev::UInt64, cost::Float64) = QueueNode( cur, prev, cost )
QueueNode(cur::Int64, prev::Int64, cost::Float64) = QueueNode( cv(cur), cv(prev), cost )


cv(x) = convert(IntGenotype,x)

function setup(n)
  #include("basins.jl")
  global ls = NKLandscape(n,2)
  global ff = lsfits(ls)
  global bcc = basincounts(basins(ls,ff),ls)
  global p1 = cv(bcc[1][1])
  global p2 = cv(bcc[2][1])
end
#setup(4)

@doc """function test_paths(n,k,num_landscapes)
Finds least-cost paths between all pairs of peaks of "num_landscapes" landscapes.
least-cost paths are computed in both directions as a check for correctness.

Results are printed for each landscape, and summarized for all landscapes.

TODO:  The "uniform cost search" algorithm is basically Djkstra's algorithm which computes the shortest (least cost)
      paths from any point to all other points.  Thus, doing separate computations for each pair is not efficient.
TODO:  This function combines testing and reporting.  These functions should be separated.
TODO:  Each landscape can be a separate thread.
"""
function test_paths(n,k,num_landscapes)
  println("fit diff weight: ",fit_diff_weight)
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
    bcc = basincounts(basins(ls,fts),ls)
    if length(bcc) > 1
      i = 1
      for i = 1:length(bcc)
        p1 = bcc[i][1]
        j = 2
        for j = (i+1):length(bcc)
          p2 = bcc[j][1]
          spw1 = shortest_path_between_peaks(p1,p2,ls,fts)
          spw2 = shortest_path_between_peaks(p2,p1,ls,fts)
          @assert length(spw1[1]) == length(spw2[1])
          @assert isapprox(spw1[2],spw2[2])
          #println("min_len: ",count_ones( p1 $ p2 ),"  length: ",length(spw1[1])-1,"  cost:",spw1[2])
          ave_length += (length(spw1[1])-1)
          ave_cost += spw1[2]
          min_length += count_ones( p1 $ p2 )  # Hamming distance between peaks plus 1
        end
      end
      denom = length(bcc)*(length(bcc)-1)/2
      ave_length /= denom
      min_length /= denom
      ave_cost /=  denom
      @printf("%4d\t%4d\t%4d\t%6.2f\t%6.2f\t%6.2f\n",n,k,length(bcc),min_length,ave_length,ave_cost)
    else
      println("single peak")
    end
    sum_ave_length += ave_length
    sum_min_length += min_length
    sum_ave_cost += ave_cost
    sumsq_ave_length += ave_length^2
    sumsq_min_length += min_length^2
    sumsq_ave_cost += ave_cost^2
  end
  std_dev_ave_length = sqrt((sumsq_ave_length - sum_ave_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_min_length = sqrt((sumsq_min_length - sum_min_length^2/num_landscapes)/(num_landscapes-1))
  std_dev_ave_cost = sqrt((sumsq_ave_cost - sum_ave_cost^2/num_landscapes)/(num_landscapes-1))
  @printf("avg_min_length:%6.2f  std_dev_min_length:%6.2f\n",sum_min_length/num_landscapes,std_dev_min_length)
  @printf("avg_ave_length:%6.2f  std_dev_ave_length:%6.2f\n",sum_ave_length/num_landscapes,std_dev_ave_length)
  @printf("avg_ave_cost:%6.2f std_dev_ave_cost:%6.2f\nn",sum_ave_cost/num_landscapes,std_dev_ave_cost)
end
    
@doc """function edge_cost(gtype1::IntGenotype, gtype2::IntGenotype, ls::NKLandscape, fits::Vector{Float64})

The cost of a path edge.  See the description at the beginning of this file.
Note the dependence of the cost on the constant "fit_diff_weight".
"""
function edge_cost(gtype1::IntGenotype, gtype2::IntGenotype, ls::NKLandscape, fits::Vector{Float64})
  fit_diff = abs(fits[gtype1+1]-fits[gtype2+1])*ls.n*fit_diff_weight
  #fit_diff = abs(fitness(gtype1,ls)-fitness(gtype2,ls))*ls.n*fit_diff_weight
  return 1.0 + fit_diff
end


@doc """function shortest_path_between_peaks(peak1::IntGenotype, peak2::IntGenotype, ls::NKLandscape, fits::Vector{Float64})

 Find the reverse least cost path between two given peaks (local fitness optima) by using A* search 
 TODO:  The "uniform cost search" algorithm is basically Djkstra's algorithm which computes the shortest (least cost)
      paths from any point to all other points.  Thus, doing separate computations for each pair is not efficient.
 TODO:  Further checking that I have the A* algorithm implemented correctly.
 TODO:  Rename variables to be more meaningful.
"""
function shortest_path_between_peaks(peak1::IntGenotype, peak2::IntGenotype, ls::NKLandscape, fits::Vector{Float64})
  closed = Set()
  queue = Base.Collections.PriorityQueue(QueueNode,Float64)
  peak1_node = QueueNode(peak1,cons(peak1,nil()),0.0)
  queue[peak1_node] =0.0
  while !isempty(queue)
    ng_node = Base.Collections.dequeue!(queue)
    if ng_node.current == peak2
      #println("final ng_node:",ng_node.current,"  ",ng_node.cost)
      len = length(ng_node.previous)
      result = zeros(IntGenotype,len)
      i = len
      for g in ng_node.previous
        result[i] = g
        i -= 1
      end
      return result,ng_node.cost
    end
    if !in(ng_node.current, closed)
      #println("ng_node:",ng_node.current,"  ",ng_node.cost)
      push!(closed, ng_node.current)
      #println("closed:",closed)
      for nbr = neighbors(ng_node.current, ls)
        #fit_diff = abs(fitness(nbr,ls)-fitness(ng_node.current,ls))*ls.n*fit_diff_weight
        #heuristic = count_ones( nbr $ peak2 )
        new_cost = edge_cost(nbr,ng_node.current,ls,fits)
        #println("cur:",ng_node.current,"  nbr:",nbr,"  cur_cost:",ng_node.cost,"   new_cost:",new_cost)
        value = get(queue,nbr,Void)
        if (value == Void) || (value < new_cost)
          nbr_node = QueueNode(nbr,cons(nbr,ng_node.previous), new_cost + ng_node.cost)
          #println("nbr:",nbr,"  nbr cost:",nbr_node.cost)
          queue[nbr_node] = new_cost + ng_node.cost
        else
          println("value >:  nbr:",nbr,"  new cost:",nbr_node.cost)
        end
      end
    end
  end
  return nil()
end

@doc """function path_cost(path::Array{IntGenotype,1},ls::NKLandscape)
  An independent computation of the fit-diff cost of a path.
  Mostly for debugging purposes.
"""
function path_cost(path::Array{IntGenotype,1},ls::NKLandscape)
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
