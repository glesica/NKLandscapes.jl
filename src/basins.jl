using DataStructures
using Base.Test

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
function basins1(ls::Landscape)
  basin_sets = IntDisjointSets(ls.a^ls.n)
  for ig::IntGenotype = 0:(ls.a^ls.n - 1)
    g = itog(ig,ls)
    fg = fitness(g,ls)
    fittest = fittest_neighbor(g,l)
    if fitness(fittest,ls) <= fg
      continue
    else
      union!(basin_sets,ig+1,gtoi(fittest,ls)+1)
    end
  end
  return basin_sets
end

@doc """basins(ls::Landscape, fits::Vector{Float64})

Returns an instance of `IntDisjointSets` that contains, as disjoint sets, the
basins of attraction (relative to a greedy adaptive walk) present in the given landscape.
Should only be applied to NKLandscapes (rather than NKpLandscapes or NKqLandscapes at this time.
"""
function basins1(ls::Landscape, fits::Vector{Float64})
  basin_sets = IntDisjointSets(ls.a^ls.n)
  for ig::IntGenotype = 0:(ls.a^ls.n - 1)
    fg = fits[ig+1]
    g = itog(ig,ls)
    #println("ig:",ig,"  ig:",ig)
    fittest = fittest_neighbor(g,ls)
    ifittest = gtoi(fittest,ls)
    if fits[ifittest+1] <= fg
      continue
    else
      union!(basin_sets,ig+1,ifittest+1)
    end
  end
  return basin_sets
end

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

# Find the fitness of the corresponding local maximum
function local_max(x,ls::Landscape)  
  gg = greedy_adaptive_walk(itog(convert(IntGenotype,x),ls),ls).history_list[:,end]
  return gtoi(gg,ls),fitness(gg,ls)
end

function basincounts(basins::DisjointSets{IntGenotype},ls::Landscape)
  c = counter(Int64)
  for i = 0:convert(UInt64,(ls.a^ls.n-1))
    #println("i:",convert(Int,i),"  fr(i):",find_root(basins,i)-1)
    push!(c,find_root(basins,i)-1)  # subtracting 1 converts to IntGenotype
  end
  basin_list = Any[]
  for k in keys(c)
    (g,fit_g) = local_max(k,ls)
    #println("k:",k,"  g:",g,"  fit_g:",fit_g)
    push!(basin_list,(g,c[k],fit_g))
  end
  #fit_local_max(x) = fitness(local_max(x),ls)
  #basin_list = map(x -> (gtoi(local_max(x),ls), c[x], fit_local_max(x)), keys(c))
  sort!(basin_list, by=x -> x[3])  # Sort by fitness
  return basin_list
end

function print_basin_summary(basin_counts)
  @printf("max_gen\tcount\tfitness\n")
  for s in basin_counts
    @printf("%6d\t%4d\t%.6f\n",s[1],s[2],s[3])
  end
end

@doc """test_basins(ls::Landscape, fits::Vector{Float64})

Tests the `basins` function by checking that the final genotype of a greedy adpative walk 
started from each possible genotype is in the same disjoint set as the starting genotype.
"""
function test_basins(ls::NKLandscape, fits::Vector{Float64})
  basins_sets = basins(ls,fits)
  local_max(x::Int64)= gtoi(greedy_adaptive_walk(itog(convert(IntGenotype,x),ls),ls).history_list[:,end],ls)
  println("local_max(1):",local_max(1))
  fit_local_max(x) = fitness(local_max(x),ls)
  for ig::IntGenotype = 0:(ls.a^ls.n - 1)
    g = itog(ig,ls)
    greedy_walk = greedy_adaptive_walk(g,ls)
    gopt = greedy_walk.history_list[:,end]
    iopt = gtoi(gopt,ls)
    @test find_root(basins_sets,ig) == find_root(basins_sets,iopt)
    print("ig:",ig,"  iopt:",iopt," fit(iopt):",fitness(iopt,ls),"  find_root(basins_sets,ig):",find_root(basins_sets,ig),"  find_root(basins_sets,iopt):",
      find_root(basins_sets,iopt))
    print("  local_max:",local_max(find_root(basins_sets,ig)))
    println("  fit_max:",fit_local_max(find_root(basins_sets,ig)))
  end
end

@doc """test_basincounts(ls::NKLandscape, fits::Vector{Float64})

Tests the `basincount` function by checking that the representative
of each basin is a local maximum.
"""
function test_basincounts(ls::NKLandscape, fits::Vector{Float64})
  basin_sets = basins(ls,fits)
  basin_counts = basincounts(basin_sets,ls)
  for b in basin_counts
    fit_b = fitness(b[1],ls)
    fn = fittest_neighbor(itog(b[1],ls),ls)
    fit_fn = fitness(fn,ls)
    @test fit_fn < fit_b
  end
end
    
      
function bc(ls::NKLandscape)
  lsf = lsfits(ls)
  lsb = basins(ls,lsf)
  return basincounts(lsb,ls)
end

