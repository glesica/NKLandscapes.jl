#= A basic genetic algorithm.
=#
using NKLandscapes

@doc """ function genetic_algorithm( ls::Landscape, pop_size::Int64, maxgens::Int64;
    selection_funct::Function=propsel!, tournsel_k::Int64=2,
    moran_sel_iters::Float64=1.0,   # Fraction of ls.n
    mutation_funct::Function=bwmutate!, mut_prob::Float64=1.0, # Mutation rate is mut_prob/ls.n
    init_genotype::Any=Void, termination_funct::Function=default_termination_funct,
    statistics_funct::Function=default_statistics_funct )

A basic generational genetic algorithm using the functions provided in NKLandscape.jl.
The non-optional arguments are the landscape, which defines the fitness function,
the populaiton size, and the maximum number of generations.
By default, the algorithm uses random initialization, proportional selection,
and bitwise mutation with rate 1/ls.n, and no crossover.  Optional arguments can
specify alternative methods of initialization, selection, and mutation, and termination.
A statistics function is called at the beginning of every generation and at termination.
"""

function genetic_algorithm( ls::Landscape, pop_size::Int64, maxgens::Int64;
    selection_funct::Function=propsel!, 
    tournsel_k::Int64=2,   # number of individuals in a tournament
    moransel_iters::Float64=1.0,   # Fraction of ls.n 
    mutation_funct::Function=bwmutate!, mut_prob::Float64=1.0, # Mutation rate is mut_prob/ls.n
    termination_funct::Function=default_termination_funct, init_genotype::Any=Void,
    statistics_funct::Function=default_statistics_funct )
  p::Population
  if typeof(init_genotype) <: Genotype
    p = constant(Population, rg, pop_size)  # Initialize to a single-genotype populaiton
  else
    p = rand(Population, ls, pop_size)
  end
  gen = 1
  while gen <= maxgens
    statistics_funct( gen, p )
    mutation_funct(p, mut_prob/ls.n)
    if selection_funct == tournsel!
      selection_funct(p, tournsel_k)
    elseif selection_funct == moransel!
      selection_funct(p, convert(Int,round(ls.n*moransel_iters)))
    else
      selection_funct(p)
    end
    if termination_funct(p) 
      break
    end
    gen += 1
  end
  statistics_funct(gen, p )
  return p
end

@doc """ function default_termination_funct( p::Population )

If this termination function is used, the algorithm runs for maxgens generations.
"""
function default_termination_funct( p::Population )
  return false
end

@doc """ function test_converged( p::Population )

Terminate if the population has converged to a single genotype.
"""
function test_converged( p::Population )
  len = length(p.genotypes)
  i = 1
  while i <= len && p.genotypes[1] == p.genotypes[i]
    i += 1
  end
  return i > len ? true : false
end

@doc """ default_statistics_funct(gen::Int64, p::Population )

Modify to collect statistics about the run.
"""
function default_statistics_funct(gen::Int64, p::Population )
  return false
end

@doc """ simple (gen::Int64, p::Population )

Report the generation number and the maximum fitness at that generation.
"""
function simple_statistics_funct(gen::Int64, p::Population )
  pfits = popfits( p )
  println("Generation: ", gen, "  Max fitness: ", maximum(pfits) )
  #println(p)
  #println(map(fitness,p.genotypes))
end

@doc """ function test_genetic_algorithm(N::Int64, K::Int64 )
  Tests the various options for function genetic_algorithm().
  Note that the use of simple_statistics_funct  does give some printed output.
"""

function test_genetic_algorithm(N::Int64, K::Int64 )
  ls = NKLandscape(N,K)
  popsize = 8
  maxgens = 10
  p0 = genetic_algorithm(ls,popsize,maxgens)
  #println("p0:",p0)

  p1 = genetic_algorithm(ls,popsize,maxgens,statistics_funct=simple_statistics_funct)
  #println("p1:",p1)
  
  p2 = genetic_algorithm(ls,popsize,maxgens,selection_funct=moransel!,moransel_iters=2.0)
  #println("p2:",p2)
  
  p3 = genetic_algorithm(ls,popsize,maxgens,selection_funct=tournsel!,tournsel_k=3)
  #println("p3:",p3)
  
  p4 = genetic_algorithm(ls,popsize,maxgens,mutation_funct=bsmutate!,mut_prob=0.5)
  #println("p4:",p4)
  
  p5 = genetic_algorithm(ls,popsize,maxgens,termination_funct=test_converged)
  #println("p5:",p5)
end

