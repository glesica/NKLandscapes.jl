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

# Random adaptive walk should terminate and move uphill

rw = random_walk(g, l)
@test rw.length > 0
for i = 2:rw.length
  @test f < fitness(rw.history[:,i], l)
end

nrw = random_walk(ng, nl)
@test nrw.length > 0
for i = 2:nrw.length
  @test nf < fitness(nrw.history[:,i], nl)
end

# Greedy adaptive walk should terminate and move uphill

gw = greedy_walk(g, l)
@test gw.length > 0
for i = 2:gw.length
  @test f < fitness(gw.history[:,i], l)
end

ngw = greedy_walk(ng, nl)
@test ngw.length > 0
for i = 2:ngw.length
  @test nf < fitness(ngw.history[:,i], nl)
end

# Reluctant adaptive walk should terminate and move uphill

rew = reluctant_walk(g, l)
@test rew.length > 0
for i = 2:rew.length
  @test f < fitness(rew.history[:,i], l)
end

nrew = reluctant_walk(ng, nl)
@test nrew.length > 0
for i = 2:nrew.length
  @test nf < fitness(nrew.history[:,i], nl)
end

