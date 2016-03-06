using Base.Collections
using DataStructures
using NKLandscapes

# Two important parameters:
const fit_diff_weight = 5.0  # Weight on fitness differences as opposed to mutational steps.
const dist_cutoff = 4    # Only paths between pairs of peaks closer to each other than this cutoff are considered

@doc """ function print_paths_summary(n,k,num_landscapes)

For each of "num_landscapes" landscapes, finds least-cost paths from each peak to all peaks of higher fitness.
Results are printed for each landscape, and summarized for all landscapes.
Note:  can use multiple processes.
"""
function print_paths_summary(n,k,num_landscapes)
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
