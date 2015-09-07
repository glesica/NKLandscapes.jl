using NK
using Base.Test

srand(0)

l = NKLandscape(10, 1)
g = random_genome(l)
f = fitness(g, l)

nl_q = NKqLandscape(10, 1, 2)
ng_q = random_genome(nl_q)
nf_q = fitness(ng_q, nl_q)

nl_p = NKpLandscape(10, 1, 0.8)
ng_p = random_genome(nl_p)
nf_p = fitness(ng_p, nl_p)

# Neighbors should differ at one locus

test_neighbors(genome, landscape) = begin
  for nbr = neighbors(genome, landscape)
    @test (genome - nbr) |> sum |> abs == 1
  end
end

test_neighbors(g, l)
test_neighbors(ng_q, nl_q)
test_neighbors(ng_p, nl_p)

# Neutral neighbors should have the same fitness

nn_q = neutral_neighbors(ng_q, nl_q)
@test size(nn_q)[1] == nl_q.n
@test size(nn_q)[2] > 0
for j = 1:size(nn_q)[2]
  @test nf_q == fitness(nn_q[:,j], nl_q)
end

nn_p = neutral_neighbors(ng_p, nl_p)
@test size(nn_p)[1] == nl_p.n
@test size(nn_p)[2] > 0
for j = 1:size(nn_p)[2]
  @test nf_p == fitness(nn_p[:,j], nl_p)
end

# Fitter neighbors should all be fitter

test_fitter_neighbors(landscape, nbrs, score) = begin
  @test size(nbrs)[1] == landscape.n
  @test size(nbrs)[2] > 0
  for j = 1:size(nbrs)[2]
    @test score < fitness(nbrs[:,j], landscape)
  end
end

fn = fitter_neighbors(g, l)
test_fitter_neighbors(l, fn, f)

nfn_q = fitter_neighbors(ng_q, nl_q)
test_fitter_neighbors(nl_q, nfn_q, nf_q)

nfn_p = fitter_neighbors(ng_p, nl_p)
test_fitter_neighbors(nl_p, nfn_p, nf_p)

# Fittest n neighbors should be all neighbors

fnn = fittest_neighbors(g, l, l.n)
@test size(fnn)[1] == l.n
@test size(fnn)[2] == number_neighbors(g, l)

nfnn_q = fittest_neighbors(ng_q, nl_q, nl_q.n)
@test size(nfnn_q)[1] == nl_q.n
@test size(nfnn_q)[2] == number_neighbors(ng_q, nl_q)

nfnn_p = fittest_neighbors(ng_p, nl_p, nl_p.n)
@test size(nfnn_p)[1] == nl_p.n
@test size(nfnn_p)[2] == number_neighbors(ng_p, nl_p)

# Fittest 1 neighbor should be fittest neighbor

fn1 = fittest_neighbor(g, l)
@test fn1 == fn[:,end]

nfn1_q = fittest_neighbor(ng_q, nl_q)
@test nfn1_q == nfn_q[:,end]

nfn1_p = fittest_neighbor(ng_p, nl_p)
@test nfn1_p == nfn_p[:,end]

# Random adaptive walk should terminate and move uphill

rw = random_walk(g, l)
@test rw.length > 0
for i = 2:rw.length
  @test f < fitness(rw.history[:,i], l)
end

nrw_q = random_walk(ng_q, nl_q)
@test nrw_q.length > 0
for i = 2:nrw_q.length
  @test nf_q < fitness(nrw_q.history[:,i], nl_q)
end

nrw_p = random_walk(ng_p, nl_p)
@test nrw_p.length > 0
for i = 2:nrw_p.length
  @test nf_p < fitness(nrw_p.history[:,i], nl_p)
end

# Greedy adaptive walk should terminate and move uphill

gw = greedy_walk(g, l)
@test gw.length > 0
for i = 2:gw.length
  @test f < fitness(gw.history[:,i], l)
end

ngw_q = greedy_walk(ng_q, nl_q)
@test ngw_q.length > 0
for i = 2:ngw_q.length
  @test nf_q < fitness(ngw_q.history[:,i], nl_q)
end

ngw_p = greedy_walk(ng_p, nl_p)
@test ngw_p.length > 0
for i = 2:ngw_p.length
  @test nf_p < fitness(ngw_p.history[:,i], nl_p)
end

# Reluctant adaptive walk should terminate and move uphill

rew = reluctant_walk(g, l)
@test rew.length > 0
for i = 2:rew.length
  @test f < fitness(rew.history[:,i], l)
end

nrew_q = reluctant_walk(ng_q, nl_q)
@test nrew_q.length > 0
for i = 2:nrew_q.length
  @test nf_q < fitness(nrew_q.history[:,i], nl_q)
end

nrew_p = reluctant_walk(ng_p, nl_p)
@test nrew_p.length > 0
for i = 2:nrew_p.length
  @test nf_p < fitness(nrew_p.history[:,i], nl_p)
end

