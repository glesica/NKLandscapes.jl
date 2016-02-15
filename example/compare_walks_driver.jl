@doc """ Driver for compare_walks_parallel.jl.
This file sets the parameters to be run.
A command line for running this program:
julia -p 8 -L "compare_walks_driver.jl" -e "run_all_jobs()"
"""
include("../src/NKLandscapes.jl")
using NKLandscapes
using Base.Test
include("compare_walks_parallel.jl")

const BOUNDARIES = ["Nearest-neighbor", "Random"]

#const N_VALUES = [60, 100, 140, 180]
const N_VALUES = [220, 240]
const K_VALUE =  8
const Q_VALUE = 0  # q for NKq landscape, set to 0 for NK landscape
const B_VALUE = 2
const RUNS = 100  # number of runs for each setting of the parameters

function run_all_jobs()
  print_header()
  trial = 1
  for n_value in N_VALUES
    jobs = [Job(trial, n_value, K_VALUE, Q_VALUE, B_VALUE) for trial = 1:RUNS ]
    run_jobs(jobs)
    trial += 1
  end
end

