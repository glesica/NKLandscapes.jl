@doc """ Compares the lengths and maximum fitnesses of random, greedy, and reluctant adaptive walks.
Outputs a table of means and standard deviations for these walks.
The parameters shown in the table are set in compare_walks_driver.jl.
See compare_walks_driver.jl for information on how to run on multiple cores.
"""
@everywhere immutable Job
  trial::Int64
  n::Int64
  k::Int64
  q::Int64  # q  is the number of integer random values for NKq landscapes (F in the paper)
            # q==0 means NK landscape, q > 0 means NKq landscape
  b::Int64
end

function run_jobs(jobs)
  results = pmap(run,jobs)
  all_results=accumulate_results(results)
  print_results(all_results)
end

@everywhere immutable JobResult
  job::Job
  rand_len::Int64
  rand_fit::Float64
  greedy_len::Int64
  greedy_fit::Float64
  reluct_len::Int64
  reluct_fit::Float64
  fiteq_len::Int64
  fiteq_fit::Float64
end

@everywhere immutable AllResults
  n::Int64
  k::Int64
  q::Int64
  rand_len::Array{Int64,1}
  rand_fit::Array{Float64,1}
  greedy_len::Array{Int64,1}
  greedy_fit::Array{Float64,1}
  reluct_len::Array{Int64,1}
  reluct_fit::Array{Float64,1}
  fiteq_len::Array{Int64,1}
  fiteq_fit::Array{Float64,1}
end

@everywhere function run(job::Job)
  if job.q == 0
    l = NKLandscape(job.n, job.k, near=job.b == 1)
  else
    l = NKqLandscape(job.n, job.k, job.q, near=job.b==1)
  end
  g = rand(Genotype, l)
  wrand = random_adaptive_walk(g, l)
  wgreedy = greedy_adaptive_walk(g, l)
  wreluct = reluctant_adaptive_walk(g, l)
  wfiteq = fitter_then_neutral_walk(g, l)
  res = JobResult(
    job,
    wrand.length,
    wrand.fitnesses[end],
    wgreedy.length,
    wgreedy.fitnesses[end],
    wreluct.length,
    wreluct.fitnesses[end],
    wfiteq.length,
    wfiteq.fitnesses[end])

  return res
end
  
@everywhere function print_header() 
  println("Length of Paths to Optima and Fitness of Optima for random, greedy, reluctant, and fitter_or_equal walks")
  println(BOUNDARIES[B_VALUE], " Interactions")
  println(RUNS, " runs per line of output")
  println("rnd:  random")
  println("grd:  greedy")
  println("rlt:  reluctant")
  println("ftq:  fitter or equal")
  println("ln:   length")
  println("ft:   fitness")
  println("N\tK\tq\tμ rndln\tσ rndln\tμ rndft\tσ rndft\tμ grdln\tσ grdln\tμ grdft\tσ grdft\tμ rltln\tσ rltln\tμ rltft\tσ rltft\tμ ftqln\tσ ftqln\tμ ftqft\tσ ftqft")
end

function accumulate_results(result_list::Array{Any,1})
  len = length(result_list)
  all_results = AllResults(
    result_list[1].job.n,
    result_list[1].job.k,
    result_list[1].job.q,
    zeros(Int64,len),
    zeros(Float64,len),
    zeros(Int64,len),
    zeros(Float64,len),
    zeros(Int64,len),
    zeros(Float64,len),
    zeros(Int64,len),
    zeros(Float64,len)
  )
  i = 1
  for res in result_list
    all_results.rand_len[i] = res.rand_len
    all_results.rand_fit[i] = res.rand_fit
    all_results.greedy_len[i] = res.greedy_len
    all_results.greedy_fit[i] = res.greedy_fit
    all_results.reluct_len[i] = res.reluct_len
    all_results.reluct_fit[i] = res.reluct_fit
    all_results.fiteq_len[i] = res.fiteq_len
    all_results.fiteq_fit[i] = res.fiteq_fit
    i += 1
  end
  return all_results
end
  
@everywhere function print_results(a::AllResults) 
    @printf("%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f", a.n, a.k, a.q, mean(a.rand_len), std(a.rand_len), mean(a.rand_fit), std(a.rand_fit))
    @printf("\t%.3f\t%.3f\t%.3f\t%.3f", mean(a.greedy_len), std(a.greedy_len), mean(a.greedy_fit), std(a.greedy_fit))
    @printf("\t%.3f\t%.3f\t%.3f\t%.3f", mean(a.reluct_len), std(a.reluct_len), mean(a.reluct_fit), std(a.reluct_fit))
    @printf("\t%.3f\t%.3f\t%.3f\t%.3f\n", mean(a.fiteq_len), std(a.fiteq_len), mean(a.fiteq_fit), std(a.fiteq_fit))
end
