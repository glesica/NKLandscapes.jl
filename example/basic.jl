include("../src/NK.jl")
using NK

srand(0)

println("NKLandscape...")
l = NKLandscape(10, 2)
println(l)

println("Random genome...")
g = rand(Genotype, l)
println(g)

println("Fitness...")
f = fitness(g, l)
println(f)

println("Neighbors...")
nbrs = all_neighbors(g, l)
for i = number_neighbors(g, l)
  nbr = nbrs[:,i]
  println(nbr)
end

println("Neutral neighbors...")
neutral_neighbors(g, l) |> println

println("Fitter neighbors...")
fitter_neighbors(g, l) |> println

println("Fittest 3 neighbors...")
fittest_neighbors(g, l, 3) |> println

println("Climbing hills...")
g0 = g
f0 = fitness(g, l)
while true
  println(g0)
  g1 = fittest_neighbor(g0, l)
  f1 = fitness(g1, l)
  if f1 > f0
    g0 = g1
    f0 = f1
  else
    break
  end
end

println("Population climbing hills...")
p = rand(Population, l, 5)
println(p)

println("Population fitnesses...")
fs = popfits(p, l)
println(fs)

println("NKp landscape...")
nl = NKpLandscape(10, 2, 0.9)
println(nl)

println("Probability neutral...")
println(prob_neutral(nl))

println("Fitness...")
nf = fitness(g, nl)
println(nf)

println("Neighbors...")
nbrs = all_neighbors(g, nl)
for i = number_neighbors(g, nl)
  nbr = nbrs[:,i]
  println(nbr)
end

println("Neutral neighbors...")
neutral_neighbors(g, nl) |> println

println("Fitter neighbors...")
fitter_neighbors(g, nl) |> println

println("Fittest 3 neighbors...")
fittest_neighbors(g, nl, 3) |> println

println("Climbing hills...")
g0 = g
f0 = fitness(g, nl)
while true
  println(g0)
  g1 = fittest_neighbor(g0, nl)
  f1 = fitness(g1, nl)
  if f1 > f0
    g0 = g1
    f0 = f1
  else
    break
  end
end

println("Population climbing hills...")
p = rand(Population, nl, 5)
println(p)

println("Population fitnesses...")
fs = popfits(p, nl)
println(fs)
