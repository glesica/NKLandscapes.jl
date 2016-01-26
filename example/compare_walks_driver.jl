include("../src/NKLandscapes.jl")
using NKLandscapes
using Base.Test
include("compare_walks_parallel.jl")

const BOUNDARIES = ["Nearest-neighbor", "Random"]

const N_VALUES = [300]
const K_VALUE = 8
const Q_VALUE = 0  # q for NKq landscape, set to 0 for NK landscape
const B_VALUE = 2
const RUNS = 100

function run_all_jobs()
  print_header()
  for n_value in N_VALUES
    jobs = [Job(trial, n_value, K_VALUE, Q_VALUE, B_VALUE) for trial = 1:RUNS ]
    run_jobs(jobs)
  end
end

