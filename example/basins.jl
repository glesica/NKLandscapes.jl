using Base.Collections
using DataStructures
using NKLandscapes


@doc """ function print_basin_summary(basin_lists)

Prints a summary of the list of basins returned by basinlists!().
Example usage:  print_basin_summary(basinlists!(basins(ls),ls))  where ls is an NKLandscape
"""
function print_basin_summary(basin_lists)
  @printf("max_gen\tcount\tfitness\n")
  for s in basin_lists
    @printf("%s\t%4d\t%.6f\n","$(bits(s.gtype.alleles)[(end - s.gtype.landscape.n + 1):end])",s.count,s.peak_fitness)
  end
end
