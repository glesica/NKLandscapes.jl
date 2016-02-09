export genotype_to_int, int_to_genotype, ifitness, fitness_array, fitness_levels_array, count_fitness_levels,
  neighbors, neutral_net, neutral_nets, list_neutral_nets, neutral_net_summary, print_neutral_nets

@doc """
Enumerates all connected neutral networks for NKqLandscapes and all
approximately neutral networks for NKLandscapes, and summarizes their
properties.  The objective is to see if there are "giant" connected neutral
networks that contain almost all genotypes at a given fitness level as
predicted by the "percolation" model of Sergey Gavrilets for uncorrelated
landscapes.  (Reference needed.) An integer representation of genotypes
is used extensively.  The functions "genotype_to_int" and
"int_to_genotype" convert between the integer representation and the
array representation used elsewhere in NKLandscapes.jl.  If the arity
of the landscape is 2, then the integer representation can be
considered as a bit-string representation, and some functions require
that the arity of the landscape is 2.

The following steps are suggested for the analysis of the connected approximately
neutral networks of an NK landscape.

```
l = NKLandscape(18,4)   # constructs the landscape
fa = fitness_array(l)   # computes all fitness and stores them in array fa (computationally expensive)
fl = fitness_levels_array(l,20,fa)  # fl is an array of fitness levels indexed on integer genotypes
FL = list_neutral_nets(l,fl)   # FL is a list of all connected neutral nets
print_neutral_nets(neutral_net_summary(FL,20))
```

"""
using NKLandscapes
using DataStructures

@doc """
Returns Genotype g converted to the equivalent base ls.a integer
"""
function genotype_to_int(g::Genotype,ls::Landscape)
  sum = g[1]-1
  for i = 2:length(g)
    sum = ls.a*sum + g[i]-1
  end
  return sum
end

@doc """
Returns base ls.a integer n converted to a Genotype (array of integers)
"""
function int_to_genotype(n::Int64,ls::Landscape)
  if ls.a == 2
    [bit == '1' ? 2 : 1 for bit in bin(n, ls.n)]
  else
    [parse(Int64,allele,ls.a)+1 for allele in base(ls.a,n,ls.n)]
  end
end

@doc """
Returns the fitness of integer genotype ig
"""
function ifitness(ig::Int64,l::Landscape)
  return fitness(int_to_genotype(ig,l),l)
end

@doc """
Returns an array of all of the fitnesses of a landscape indexed by integer genotypes.
This function is computationally expensive.
TODO:  Figure out how to parallelize this function.  Note that doing a pmap over genotypes
  seems to be too fine-grained.  Maybe break the list of genotypes into blocks of 1000?
"""
function fitness_array(l::Landscape)
  fitness_array = zeros(Float64,l.a^l.n)
  for i = 0:l.a^l.n-1
    fitness_array[i+1] = ifitness(i,l)
  end
  return fitness_array
end

@doc """
Returns fitness_levels which is an array of the fitness level of each integer genotype of the landscape
Assumes that the interval [0,1] is subdivided into  num_intervals  subintervals.
fitness_levels[ig] is the index of the interval for integer genotype ig
The fitness array may be given as the third argument
"""
function fitness_levels_array(l::Landscape,num_intervals::Int, fa::Array{Float64,1}=Float64[])
  fitness_levels = zeros(Int64,l.a^l.n)
  if typeof(l)==NKLandscapes.NKqLandscape
    eps_inc = eps()  # slightly reduce the lower bound on intervals for NKq landscapes
  else
    eps_inc = 0.0
  end
  if fa == Float64[]
    for i = 0:l.a^l.n-1
      fitness_levels[i+1] = floor(Int,(ifitness(i,l)+eps_inc)*num_intervals)
    end
  else
    for i = 0:l.a^l.n-1
      fitness_levels[i+1] = floor(Int,(fa[i+1]+eps_inc)*num_intervals)
    end
  end
  return fitness_levels
end

@doc """
Returns an array of the number of genotypes in each fitness level
"""
function count_fitness_levels(l::Landscape,num_intervals::Int,fl::Array{Int,1}=Int[])
  if fl == Int[]
    fl = fitness_levels_array(l,num_intervals)
  end
  counts = zeros(Int64,num_intervals)
  for i = 0:(l.a^l.n-1)
    counts[fl[i+1]+1] += 1
  end
  return counts
end

@doc """
Returns an array of the of integer genotype neighbors of integer genotype ig
Assumes landscape l has arity 2, i. e.,  l.a==2.
"""
function neighbors( ig::Int, l::Landscape )
  @assert l.a == 2
  nbrs_list = zeros(Int,l.n)
  single_bit = convert(UInt64,0x1)
  for i = 0:l.n-1
    nbrs_list[i+1] = ig $ single_bit
    single_bit <<= 1
  end
  nbrs_list
end

@doc """
Returns the connected neutral network of integer genotype ig relative to the given fitness_levels array.
In other words, two integer genotypes are neutral neighbors if they have the same fitness level.
"""
function neutral_net(ig::Int,l::Landscape,fitness_levels::Array{Int,1})
  @assert l.a == 2
  ig0 = ig
  f0 = fitness_levels[ig+1]
  closed = Set()
  stack = Stack(Int)
  push!(stack,ig0)
  while !isempty(stack)
    ig = pop!(stack)
    if !in(ig,closed)
      push!(closed,ig)
      nbrs = filter(n->fitness_levels[n+1]==f0,neighbors(ig,l))
      for n in nbrs
        push!(stack,n)
      end
    end
  end
  return collect(closed)
end
  
@doc """
Adds the Ints in int_list to the IntDisjointSets data structure S
"""
function add_ints_to_disjoint_sets!(S::IntDisjointSets,int_list::Array{Any,1})
  for i = 2:length(int_list)
    union!(S,int_list[i-1]+1,int_list[i]+1)
  end
  S
end

@doc """
For this and succeeding functions, a "neutral net" is a path connected set of genotypes with the same
    fitness level.  Two genotypes are "connected" if they are neighbors (i. e., differ at one locus).
Returns the IntDisjointSets representation of all of the connected neutral nets of landscape l.
"""
function neutral_nets(l::Landscape,fitness_levels::Array{Int,1})
  @assert l.a == 2
  S = IntDisjointSets(l.a^l.n)
  for i = 0:(l.a^l.n-1)
    if (find_root(S,i+1) != i+1) || (S.ranks[i+1] > 0)
      continue   # i is already in S
    end
    nbrs = neutral_net(i,l,fitness_levels)
    if length(nbrs) > 0
      push!(nbrs,i)  # add i to nbrs
    else
      nbrs = Int64[i]   # Array containing only i
    end
    add_ints_to_disjoint_sets!(S,nbrs)
  end
  return S
end

@doc """

Returns a list of all of the neutral nets (relative to fitness_levels) of landscape l.
The returned list is sorted by the size (number of genotypes) of the neutral net.
Each element of the list is a triple with the following components:
*  An integer genotype representative of the neutral net
*  The number of genotypes in the neutral net
*  The integer fitness level of the neutral net
The returned list is sorted on the number of genotypes in the neutral nets.
"""
function list_neutral_nets(l::Landscape,fitness_levels::Array{Int,1})
  @assert l.a == 2
  C = counter(Int)
  S = neutral_nets(l,fitness_levels)
  for i = 0:(l.a^l.n-1)
    push!(C,find_root(S,i+1)-1)
  end
  nn_list = map(x->(x,C[x],fitness_levels[x+1]), keys(C))
  sort!(nn_list,lt=(x,y)->x[2]<y[2])   # sort by size of neutral net
  return nn_list
end

@doc """
Returns a list of information about the connected neutral nets at each non-empty fitness level
Each list element is a 5-tuple with the following components 
*  The fit level (empty fit levels are not included)
*  The lower bound of the fit level
*  The number of genotypes at the fit level
*  The maximum size of a neutral net at the fit level
*  The ratio of maximum size neutral net over the number of genotypes
Note that if this ratio is close to 1.0, then a high proportion of the genotypes at the fitness level 
  are in a "giant component" neutral net.
"""
function neutral_net_summary(nn_list,num_intervals)
  fit_increment = 1.0/num_intervals
  nn3 = map(x->x[3],nn_list)
  lb = minimum(nn3)
  ub = maximum(nn3)
  summary = []
  for i = lb:ub
    filtered_nn = [y[2] for y in filter(x->x[3]==i,nn_list)]
    if length(filtered_nn) > 0
      max_nn = maximum(filtered_nn)
      sum_nn = sum(filtered_nn)
      push!(summary,(i, i*fit_increment, sum_nn, max_nn, max_nn/sum_nn))
    end
  end
  return summary
end

@doc """
Prints the the list of fitness levels returned by neutral_net_summary in CSV format.
"""
function print_neutral_nets(summary)
  @printf("fit lev\tfit_lb\t total count\t   max count\tpct max\n")
  for s in summary
    @printf("%2d\t%.3f\t%12d\t%12d\t%.5f\n",s[1],s[2],convert(Int,s[3]),convert(Int,s[4]),s[5])
  end
end
