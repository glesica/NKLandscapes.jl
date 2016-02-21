@doc """
Enumerates all connected neutral networks for NKqLandscapes and all
approximately neutral networks for NKLandscapes, and summarizes their
properties.  The objective is to see if there are "giant" connected neutral
networks that contain almost all genotypes at a given fitness level as
predicted by the "percolation" model of Sergey Gavrilets for uncorrelated
landscapes.  (Reference needed.) An integer representation of genotypes
is used extensively.  The functions `gtoi` and
`itog` convert between the integer representation and the
array representation used elsewhere in NKLandscapes.jl.  If the arity
of the landscape is 2, then the integer representation can be
considered as a bit-string representation, and some functions require
that the arity of the landscape is 2.
"""
using NKLandscapes
using DataStructures

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
function neutral_net_summary(counts, num_intervals)
  fit_increment = 1.0/num_intervals
  nn3 = map(x->x[3],counts)
  lb = minimum(nn3)
  ub = maximum(nn3)
  summary = []
  for i = lb:ub
    filtered_nn = [y[2] for y in filter(x->x[3]==i,counts)]
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
