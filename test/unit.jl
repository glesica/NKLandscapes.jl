using NK
using FactCheck

srand(1)

facts("Landscapes") do
  l = NKLandscape(10, 1)
  g = random_genome(l)
  f = fitness(g, l)

  nl_q = NKqLandscape(10, 1, 2)
  ng_q = random_genome(nl_q)
  nf_q = fitness(ng_q, nl_q)

  nl_p = NKpLandscape(10, 1, 0.8)
  ng_p = random_genome(nl_p)
  nf_p = fitness(ng_p, nl_p)

  context("Neighbors should differ at one locus") do
    function test_neighbors(genome, landscape)
      nbrs = all_neighbors(genome, landscape)
      for i = 1:number_neighbors(genome, landscape)
        nbr = nbrs[:, i]
        @fact (genome - nbr) |> sum |> abs --> 1
      end
    end

    test_neighbors(g, l)
    test_neighbors(ng_q, nl_q)
    test_neighbors(ng_p, nl_p)
  end

  context("Neutral neighbors should have the same fitness") do
    function test_neutral_neighbors(genome, landscape)
      nbrs = neutral_neighbors(genome, landscape)
      score = fitness(genome, landscape)
      @fact size(nbrs)[1] --> landscape.n
      @fact size(nbrs)[2] --> greater_than(0)
      for j = 1:size(nbrs)[2]
        @fact fitness(nbrs[:,j], landscape) --> score
      end
    end

    test_neutral_neighbors(ng_q, nl_q)
    test_neutral_neighbors(ng_p, nl_p)
  end

  context("Fitter neighbors should all be fitter") do
    function test_fitter_neighbors(nbrs, landscape, score)
      @fact size(nbrs)[1] --> landscape.n
      @fact size(nbrs)[2] --> greater_than(0)
      for j = 1:size(nbrs)[2]
        @fact fitness(nbrs[:,j], landscape) --> greater_than(score)
      end
    end

    fn = fitter_neighbors(g, l)
    test_fitter_neighbors(fn, l, f)

    nfn_q = fitter_neighbors(ng_q, nl_q)
    test_fitter_neighbors(nfn_q, nl_q, nf_q)

    nfn_p = fitter_neighbors(ng_p, nl_p)
    test_fitter_neighbors(nfn_p, nl_p, nf_p)
  end

  context("Fittest n neighbors should be all neighbors") do
    fnn = fittest_neighbors(g, l, l.n)
    @fact size(fnn)[1] --> l.n
    @fact size(fnn)[2] --> number_neighbors(g, l)

    nfnn_q = fittest_neighbors(ng_q, nl_q, nl_q.n)
    @fact size(nfnn_q)[1] --> nl_q.n
    @fact size(nfnn_q)[2] --> number_neighbors(ng_q, nl_q)

    nfnn_p = fittest_neighbors(ng_p, nl_p, nl_p.n)
    @fact size(nfnn_p)[1] --> nl_p.n
    @fact size(nfnn_p)[2] --> number_neighbors(ng_p, nl_p)
  end

  context("Fittest 1 neighbor should be the fittest neighbor") do
    function test_fittest_neighbor(genome, landscape)
      nbrs = fitter_neighbors(genome, landscape)
      nbr1 = fittest_neighbor(genome, landscape)
      @fact nbr1 --> nbrs[:,end]
    end

    test_fittest_neighbor(g, l)
    test_fittest_neighbor(ng_q, nl_q)
    test_fittest_neighbor(ng_p, nl_p)
  end

  context("Adaptive walks") do
    function test_adaptive_walk(walk_function, genome, landscape)
      walk = walk_function(genome, landscape)
      @fact walk.length --> greater_than(0)
      for i = 2:walk.length
        @fact fitness(walk.history[:,i], landscape) --> greater_than(fitness(genome, landscape))
      end
    end

    context("Random adaptive walk should terminate and move uphill") do
      test_adaptive_walk(random_walk, g, l)
      test_adaptive_walk(random_walk, ng_q, nl_q)
      test_adaptive_walk(random_walk, ng_p, nl_p)
    end

    context("Greedy adaptive walk should terminate and move uphill") do
      test_adaptive_walk(greedy_walk, g, l)
      test_adaptive_walk(greedy_walk, ng_q, nl_q)
      test_adaptive_walk(greedy_walk, ng_p, nl_p)
    end

    context("Reluctant adaptive walk should terminate and move uphill") do
      test_adaptive_walk(reluctant_walk, g, l)
      test_adaptive_walk(reluctant_walk, ng_q, nl_q)
      test_adaptive_walk(reluctant_walk, ng_p, nl_p)
    end
  end
end

