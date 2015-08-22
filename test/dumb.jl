using NK
using Base.Test

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

println("Random neighbor...")
for n = random_neighbor(g, l)
  println(n)
end

println("Neutral neighbors...")
for n = neutral_neighbors(g, l)
  println(n)
end

println("Better neighbors...")
for n = better_neighbors(g, l)
  println(n)
end

println("Best 3 neighbors...")
for n = best_neighbors(g, l, 3)
  println(n)
end

println("Climbing hills...")
g0 = g
f0 = fitness(g, l)
while true
  println(g0)
  g1 = best_neighbor(g0, l)
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

println("Fitness...")
nf = fitness(g, nl)
println(nf)

println("Neighbors...")
for n = neighbors(g, nl)
  println(n)
end

println("Random neighbor...")
for n = random_neighbor(g, nl)
  println(n)
end

println("Neutral neighbors...")
for n = neutral_neighbors(g, nl)
  println(n)
end

println("Better neighbors...")
for n = better_neighbors(g, nl)
  println(n)
end

println("Best 3 neighbors...")
for n = best_neighbors(g, nl, 3)
  println(n)
end

println("Climbing hills...")
g0 = g
f0 = fitness(g, nl)
while true
  println(g0)
  g1 = best_neighbor(g0, nl)
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