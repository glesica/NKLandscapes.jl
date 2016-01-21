# Investigate the properties of the connected neutral networks of an NKqLandscape
include("../src/NKLandscapes.jl")
using NKLandscapes
using DataStructures

# Returns Genotype g converted to the equivalent base ls.a integer
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

# Returns the fitness of integer genotype ig
function ifitness(ig::Int64,l::Landscape)
  return fitness(int_to_genotype(ig,l),l)
end

# Returns the sequence of integer genotypes of a neutral walk starting at genotype g
function neutral_walk_to_ints(g::Genotype, l::Landscape )
  w = neutral_walk(g,l)
  return map(x->genotype_to_int(x,l), [w.history_list[:,i] for i = 1:size(w.fitnesses)[1]])
end

# Returns the integer genotypes of the neutral neighbors of genotype g
function neutral_neighbors_to_ints(g::Genotype, l::Landscape )
  nbrs = neutral_neighbors(g,l)
  return map(x->genotype_to_int(x,l), [nbrs[:,i] for i = 1:size(nbrs)[2]])
end

function add_ints_to_disjoint_sets!(S::IntDisjointSets,int_list::Array{Any,1})
  for i = 2:length(int_list)
    union!(S,int_list[i-1]+1,int_list[i]+1)
  end
  S
end

# Returns an IntDisjointSets representation of the neutral nets of the landscape
function build_neutral_nets(l::Landscape)
  S = IntDisjointSets(l.a^l.n)
  for i = 0:(l.a^l.n-1)
    if (find_root(S,i+1) != i+1) || (S.ranks[i+1] > 0)
      continue   # i is already in S
    end
    nbrs= neutral_net(i,l)
    if length(nbrs) > 0
      push!(nbrs,i)  # add i to nbrs
    else
      nbrs = Int64[i]
    end
    add_ints_to_disjoint_sets!(S,nbrs)
  end
  return S
end

# Returns a list of the integer genotypes in the connected neutral net of integer genotype ig
# Does a depath-first traversal of the neutral net
function neutral_net(ig::Int,l::Landscape)
  if  ig >= l.a^l.n
   error("integer genotype is too large for landscape") 
  end
  if ig < 0
   error("integer genotype is negative") 
  end
  ig0 = ig
  f0 = ifitness(ig0,l)
  closed = Set()
  stack = Stack(Int)
  push!(stack,ig0)
  while !isempty(stack)
    ig = pop!(stack)
    if !in(ig,closed)
      push!(closed,ig)
      nbrs = neutral_neighbors_to_ints(int_to_genotype(ig,l),l)
      for n in nbrs
        push!(stack,n)
      end
    end
  end
  return collect(closed)
end

# Returns a list of all connected neutral networks of the landscape.
# Each neutral net is represented by a triple whose elements are:
#   1)  an integer genotype representative of the network
#   2)  the number of genotypes in the neutral network
#   3)  the fitness of the neutral network
function list_neutral_nets(l::Landscape)
  C = counter(Int)
  S = build_neutral_nets(l)
  for i = 0:(l.a^l.n-1)
    push!(C,find_root(S,i+1)-1)
  end
  return map(x->(x,C[x],ifitness(x,l)), keys(C))
end
