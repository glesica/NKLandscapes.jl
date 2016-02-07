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
fl = fitness_levels_array(l,fa,20)  # fl is an array of fitness levels indexed on integer genotypes
FL = list_neutral_nets(l,fl)   # FL is a list of all connected neutral nets
print_summary(summary_list(FL,20))
```

TODO: interpret the output of "summarize".  
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

# Returns base ls.a integer n converted to a Genotype (array of integers)
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

function add_ints_to_disjoint_sets!(S::IntDisjointSets,int_list::Array{Any,1})
  for i = 2:length(int_list)
    union!(S,int_list[i-1]+1,int_list[i]+1)
  end
  S
end

function fitness_array(l::Landscape)
  fitness_array = zeros(Float64,l.a^l.n)
  for i = 0:l.a^l.n-1
    fitness_array[i+1] = ifitness(i,l)
  end
  return fitness_array
end

# fitness_levels is an array of the fitness level of each integer genotype of the landscape
# Assumes that the interval [0,1] is subdivided into  num_intervals  subintervals.
# fitness_levels[ig] is the index of the interval for integer genotype ig
function fitness_levels_array(l::Landscape,num_intervals::Int)
  fitness_levels = zeros(Int64,l.a^l.n)
  if typeof(l)==NKLandscapes.NKqLandscape
    eps_inc = eps()
  else
    eps_inc = 0.0
  end
  for i = 0:l.a^l.n-1
    fitness_levels[i+1] = floor(Int,(ifitness(i,l)+eps_inc)*num_intervals)
  end
  return fitness_levels
end

# fitness_levels is an array of the fitness level of each integer genotype of the landscape
# Assumes that the interval [0,1] is subdivided into  num_intervals  subintervals.
# fitness_levels[ig] is the index of the interval for integer genotype ig
function fitness_levels_array(l::Landscape,fa::Array{Float64,1},num_intervals::Int)
  fitness_levels = zeros(Int64,l.a^l.n)
  if typeof(l)==NKLandscapes.NKqLandscape
    eps_inc = eps()
  else
    eps_inc = 0.0
  end
  for i = 0:l.a^l.n-1
    fitness_levels[i+1] = floor(Int,(fa[i+1]+eps_inc)*num_intervals)
  end
  return fitness_levels
end

function count_fitness_levels(l::Landscape,num_intervals::Int)
  fit_levels = fitness_levels_array(l,num_intervals)
  counts = zeros(Int64,num_intervals)
  for i = 0:l.a^l.n-1
    counts[fit_levels[i+1]] += 1
  end
  return counts
end

function count_fitness_levels(l::Landscape,fl::Array{Any,1},num_intervals::Int)
  counts = zeros(Int64,num_intervals)
  for i = 0:l.a^l.n-1
    counts[fl[i+1]] += 1
  end
  return counts
end

# assumes landscape l has arity 2, i. e.,  l.a==2.
function neighbors( ig::Int, l::Landscape )
  @assert l.a == 2
  nbrs_list = zeros(Int,l.n)
  single_bit = convert(UInt64,0x1)
  for i = 0:l.n-1
    #@printf("%#06X\n",single_bit)
    nbrs_list[i+1] = ig $ single_bit
    single_bit <<= 1
  end
  nbrs_list
end

# fitness_levels is an array of the fitness level of each integer genotype of the landscape
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

function list_neutral_nets(l::Landscape,fitness_levels::Array{Int,1})
  @assert l.a == 2
  C = counter(Int)
  S = neutral_nets(l,fitness_levels)
  for i = 0:(l.a^l.n-1)
    #if fitness_levels[i+1] >= lower_bound
      push!(C,find_root(S,i+1)-1)
    #end
  end
  nn_list = map(x->(x,C[x],fitness_levels[x+1]), keys(C))
  sort!(nn_list,lt=(x,y)->x[2]<y[2])   # sort by size of neutral net
  return nn_list
end

function summary_list(nn_list,num_intervals)
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

function print_summary(summary)
  @printf("fit lev\tfit_lb\t total count\t   max count\tpct max\n")
  for s in summary
    @printf("%2d\t%.3f\t%12d\t%12d\t%.5f\n",s[1],s[2],convert(Int,s[3]),convert(Int,s[4]),s[5])
  end
end
