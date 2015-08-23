include("../src/NK.jl")
using NK

srand(0)

println("Landscape...")
l = Landscape(10, 2)
println(l)

println("Random genome...")
g = random_genome(l)
println(g)

println("Fitness...")
f = fitness(g, l)
println(f)

println("Neighbors...")
for n = neighbors(g, l)
  println(n)
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
p = random_population(5, l)
println(p)

println("Population fitnesses...")
fs = fitnesses(p, l)
println(fs)

println("Neutral landscape...")
nl = Landscape(10, 2, 0.9)
println(nl)

println("Probability neutral...")
println(prob_neutral(nl))

println("Fitness...")
nf = fitness(g, nl)
println(nf)

println("Neighbors...")
for n = neighbors(g, nl)
  println(n)
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
p = random_population(5, nl)
println(p)

println("Population fitnesses...")
fs = fitnesses(p, nl)
println(fs)
