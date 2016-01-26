# Investigate the properties of the connected neutral networks of an NKqLandscape
include("../src/NKLandscapes.jl")
using NKLandscapes
using DataStructures

@doc """ Returns Genotype g converted to the equivalent base ls.a integer
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

@doc """ Returns the fitness of integer genotype ig
"""
function ifitness(ig::Int64,l::Landscape)
  return fitness(int_to_genotype(ig,l),l)
end

@doc """ Returns the sequence of integer genotypes of a neutral walk starting at genotype g
"""
function neutral_walk_to_ints(g::Genotype, l::Landscape )
  w = neutral_walk(g,l)
  return map(x->genotype_to_int(x,l), [w.history_list[:,i] for i = 1:size(w.fitnesses)[1]])
end

@doc """ Returns the integer genotypes of the neutral neighbors of genotype g
"""
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

@doc """ Returns an IntDisjointSets representation of the neutral nets of the landscape
"""
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

@doc """ Returns a list of the integer genotypes in the connected neutral net of integer genotype ig
Does a depath-first traversal of the neutral net
"""
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

@doc """ Returns a triple where the first element nn_list of the triple is a list of all connected neutral 
   networks of the landscape.
Each neutral net of nn_list is represented by a triple whose elements are:
   1)  an integer genotype representative of the network
   2)  the number of genotypes in the neutral network
   3)  the fitness of the neutral network
The second element S of the returned triple is the disjoint sets representation of the collection of neutral nets.
The third element C of the returned triple is a counter which contains the count of each root of S
"""
function list_neutral_nets(l::Landscape)
  C = counter(Int)
  S = build_neutral_nets(l)
  for i = 0:(l.a^l.n-1)
    push!(C,find_root(S,i+1)-1)
  end
  nn_list = map(x->(x,C[x],ifitness(x,l)), keys(C))
  sort!(nn_list,lt=(x,y)->x[3]<y[3])   # sort by fitness
  return nn_list,S,C
end

@doc""" Writes two CSV files that show all neutral nets.
Onc CSV file is sorted by the size of the neutral nets.
The other CSV file is sorted by the fitness of the neutral nets.
nn_lst is the list of neutral nets that is computed by the list_neutral_nets() function
"""
function write_neutral_nets(l::Landscape,nn_lst=Void)
  if nn_lst == Void
    (nn_lst,nn_d_sets,nn_counter)  = list_neutral_nets(l)  # nn_d_sets is the disjoint sets representation of the neutral nets
  end
  lst_count_sort = sort(nn_lst,lt=(x,y)->x[2]<y[2])
  fcs = open("list_nets/list_net_N"*string(l.n)*"_K"*string(l.k)*"_q"*string(l.q)*"_cnt_"*string(Dates.today())*".csv",
    false,true,true,false,false)
  for i = 1:length(nn_lst)
    println(fcs,lst_count_sort[i][2],",",lst_count_sort[i][3])
  end
  close(fcs)
  lst_fit_sort = sort(nn_lst,lt=(x,y)->x[3]<y[3])
  ffs = open("list_nets/list_net_N"*string(l.n)*"_K"*string(l.k)*"_q"*string(l.q)*"_fit_"*string(Dates.today())*".csv",
    false,true,true,false,false)
  for i = 1:length(nn_lst)
    println(ffs,lst_fit_sort[i][2],",",lst_fit_sort[i][3])
  end
  close(ffs)
end
