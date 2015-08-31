using NK
using Base.Test

srand(0)

l = Landscape(5, 1)
g = random_genome(l)
f = fitness(g, l)

nl = Landscape(5, 1, 0.8)
ng = random_genome(nl)
nf = fitness(ng, nl)

# Neighbors should differ at one locus

for nbr = neighbors(g, l)
  @test (g - nbr) |> sum |> abs == 1
end

for nbr = neighbors(ng, nl)
  @test (ng - nbr) |> sum |> abs == 1
end

# Neutral neighbors should have the same fitness

nn = neutral_neighbors(ng, nl)
@test size(nn)[1] == nl.n
@test size(nn)[2] > 0
for j = 1:size(nn)[2]
  @test nf == fitness(nn[:,j], nl)
end

# Fitter neighbors should all be fitter

fn = fitter_neighbors(g, l)
@test size(fn)[1] == l.n
@test size(fn)[2] > 0
for j = 1:size(fn)[2]
  @test f < fitness(fn[:,j], l)
end

nfn = fitter_neighbors(ng, nl)
@test size(nfn)[1] == nl.n
@test size(fn)[2] > 0
for j = 1:size(nfn)[2]
  @test nf < fitness(nfn[:,j], l)
end

# Fittest n neighbors should be all neighbors

fnn = fittest_neighbors(g, l, l.n)
@test size(fnn)[1] == l.n
@test size(fnn)[2] == number_neighbors(g, l)

nfnn = fittest_neighbors(ng, nl, nl.n)
@test size(nfnn)[1] == nl.n
@test size(nfnn)[2] == number_neighbors(ng, nl)

# Fittest 1 neighbor should be fittest neighbor

fn1 = fittest_neighbor(g, l)
@test fn1 == fn[:,end]

nfn1 = fittest_neighbor(ng, nl)
@test nfn1 == nfn[:,end]

