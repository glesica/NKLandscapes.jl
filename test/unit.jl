using NK
using Base.Test

srand(0)

l = Landscape(10, 2)
g = random_genome(l)
f = fitness(g, l)

nl = Landscape(10, 2, 0.9)
ng = random_genome(nl)
nf = fitness(g, nl)

# Neighbors should differ at one locus

for nbr = neighbors(g, l)
  @test (g - nbr) |> sum |> abs == 1
end

for nbr = neighbors(ng, nl)
  @test (ng - nbr) |> sum |> abs == 1
end

# Neutral neighbors should have the same fitness

nn = neutral_neighbors(ng, l)
for j = 1:size(nn)[2]
  @test nf == fitness(nn[:,j], nl)
end

# Fitter neighbors should all be fitter

fn = fitter_neighbors(g, l)
for j = 1:size(fn)[2]
  @test f < fitness(fn[:,j], l)
end

nfn = fitter_neighbors(ng, nl)
for j = 1:size(nfn)[2]
  @test nf < fitness(nfn[:,j], l)
end

# Fittest n neighbors should be all neighbors

fnn = fittest_neighbors(g, l, l.n)
@test size(fnn)[2] == number_neighbors(g, l)

nfnn = fittest_neighbors(ng, nl, nl.n)
@test size(nfnn)[2] == number_neighbors(ng, nl)

