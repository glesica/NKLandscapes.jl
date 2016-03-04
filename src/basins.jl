export Basin, basins, basincounts, print_basin_summary
using DataStructures
using Base.Test

# Needs considerable "clean up", such as deleting debugging statements, improving documentation, deleting debugging functions, etc.

type Basin
  gtype::Genotype
  count::Int64
  peak_fitness::Float64
end

function ranks(s::DisjointSets{IntGenotype})
  return s.internal.ranks
end

function parents(s::DisjointSets{IntGenotype})
  return s.internal.parents
end

@doc """basins(ls::Landscape)

Returns an instance of `IntDisjointSets` that contains, as disjoint sets, the
basins of attraction (relative to a greedy adaptive walk) present in the given landscape.
Should only be applied to NKLandscapes (rather than NKpLandscapes or NKqLandscapes at this time.
"""
function basins(ls::Landscape)
  basin_sets = IntDisjointSets(ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    fg = fitness(g)
    fittest = fittest_neighbor(g)
    #println("g:",g,"  fg:",fg,"  fittest:",fittest)
    if fitness(fittest) <= fg
      continue
    else
      union!(basin_sets,g.alleles+1,fittest.alleles+1)
    end
  end
  return basin_sets
end

function basins(ls::Landscape, fits::Vector{Float64})
  basin_sets = IntDisjointSets(ls.a^ls.n)
  for i = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    if fits != Void
      fg = fits[i+1]
    else
      fg = fitness(g)
    end
    fittest = fittest_neighbor(g)
    #println("g:",g,"  fg:",fg,"  fittest:",fittest)
    if length(fits) > 0
      f_fittest = fits[fittest.alleles+1]
    else
      f_fittest = fitness(g)
    end
    if f_fittest <= fg
      continue
    else
      union!(basin_sets,g.alleles+1,fittest.alleles+1)
    end
  end
  return basin_sets
end

#=
function basins(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = DisjointSets{IntGenotype}(collect(0:convert(IntGenotype,(ls.a^ls.n-1))))
  for ig::IntGenotype = 0:(ls.a^ls.n - 1)
    fg = fits[ig+1]
    g = itog(ig,ls)
    fittest = fittest_neighbor(g,ls)
    ifittest = gtoi(fittest,ls)
    #print("ig:",ig,"  fg:",fg,"  ifittest:",ifittest,"  fit(ifittest):",fits[ifittest+1])
    if fits[ifittest+1] <= fg
      #println(" no union")
      continue
    else
      union!(basin_sets,ig,ifittest)
      #println("    union")
    end
  end
  return basin_sets
end
=#

# Find the fitness of the corresponding local maximum
function local_max(g::Genotype)  
  return greedy_adaptive_walk(g).history_list[end]
end

@doc """function change_root!(s::IntDisjointSets, new_root::Int64)

Changes the root of the set containing "new_root" to be "new_root". 
Does not change the membership of any set.
d.ranks[new_root] may be too large, but this won't change the correctness of any subsequent operation.
"""
function change_root!(s::IntDisjointSets, new_root::UInt128)
  new_root = convert(Int64,new_root)
  old_root::Int64 = find_root(s,new_root)
  s.parents[old_root] = new_root
  s.parents[new_root] = new_root
  #println("new root:",new_root,"  ",s.ranks[new_root],"  old_root:",old_root,"  ",s.ranks[old_root])
  s.ranks[new_root] = s.ranks[old_root]+1
  println("new root:",new_root,"  ",s.ranks[new_root],"  old_root:",old_root,"  ",s.ranks[old_root])
  println("find_root(new_root):",find_root(s,new_root))
end

@doc """function basincounts(basins::DisjointSets{IntGenotype},ls::Landscape)

"""

function basincounts!(basins::IntDisjointSets,ls::Landscape)
  c = counter(Int64)
  for i = 0:(ls.a^ls.n-1)
    #println("i:",i,"  fr(i):",find_root(basins,i+1)-1)
    push!(c,find_root(basins,i+1)-1)  
  end
  basin_list = Basin[]
  for k in keys(c)
    g = local_max(Genotype(k,ls))
    fit_g = fitness(g)
    change_root!(basins,g.alleles+1)
    println("k:",k,"  g:",g,"  fit_g:",fit_g,"  g.alleles:",g.alleles)
    println("frg:",find_root(basins,g.alleles+1),"  fit:",fitness(Genotype(find_root(basins,g.alleles+1),ls)))
    push!(basin_list,Basin(g,c[k],fit_g))
  end
  sort!(basin_list, by=x -> x.peak_fitness)  # Sort by fitness
  return basin_list
end

function print_basin_summary(basin_counts)
  @printf("max_gen\tcount\tfitness\n")
  for s in basin_counts
    @printf("%s\t%4d\t%.6f\n","$(bits(s.gtype.alleles)[(end - s.gtype.landscape.n + 1):end])",s.count,s.peak_fitness)
  end
end

@doc """test_basins(ls::Landscape, fits::Vector{Float64})

Tests the `basins` function by checking that the final genotype of a greedy adpative walk 
started from each possible genotype is in the same disjoint set as the starting genotype.
TODO:  Move this function to test/basins.jl
"""
function test_basins(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = basins(ls,fits)
  basin_counts = basincounts!(basin_sets,ls)
  local_max(x::Int64)= greedy_adaptive_walk(Genotype(x,ls)).history_list[end]
  fit_local_max(x) = fitness(local_max(x))
  for i::Int64 = 0:(ls.a^ls.n - 1)
    g = Genotype(i,ls)
    greedy_walk = greedy_adaptive_walk(g)
    gopt = greedy_walk.history_list[end]
    @test find_root(basin_sets,i+1) == find_root(basin_sets,gopt.alleles+1)
    @test find_root(basin_sets,i+1) == convert(Int64,gopt.alleles+1)
    #print("ig:",ig,"  iopt:",iopt," fit(iopt):",fitness(iopt,ls),"  find_root(basin_sets,ig):",find_root(basin_sets,ig),"  find_root(basin_sets,iopt):",
    #  find_root(basin_sets,iopt))
    #print("  local_max:",local_max(find_root(basin_sets,i+1)))
    #println("  fit_max:",local_max(find_root(basin_sets,g.alleles+1)))
  end
  # Chekc that the representative of each basin is a local maximum.
  for b in basin_counts
    fit_b = fitness(b.gtype)
    fn = fittest_neighbor(b.gtype)
    fit_fn = fitness(fn)
    @test fit_fn < fit_b
  end
end

@doc """test_basincounts!(ls::NKLandscape, fits::Vector{Float64})

Tests the `basincount` function by checking that the representative
of each basin is a local maximum.
TODO:  Move this function to test/basins.jl
"""
function test_basincounts!(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = basins(ls,fits)
  basin_counts = basincounts!(basin_sets,ls)
  for b in basin_counts
    fit_b = fitness(b.gtype)
    fn = fittest_neighbor(b.gtype)
    fit_fn = fitness(fn)
    @test fit_fn < fit_b
  end
end
    
function test_basincounts!(ls::NKLandscape)
  basin_sets = basins(ls)
  basin_counts = basincounts!(basin_sets,ls)
  for b in basin_counts
    fit_b = fitness(b.gtype)
    fn = fittest_neighbor(b.gtype)
    fit_fn = fitness(fn)
    @test fit_fn < fit_b
  end
end
    
# temporary utility function for debugging.  Should be eventually deleted.      
function bc(ls::NKLandscape)
  lsf = lsfits(ls)
  lsb = basins(ls,lsf)
  return basincounts!(lsb,ls)
end

