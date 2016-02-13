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
