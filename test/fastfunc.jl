import NKLandscapes
const NK = NKLandscapes
using FactCheck

n = 3
k = 1
a = 2
links = [1 2 3;
         2 3 1]
contribs = NK.Contribs(
  # Locus 1
  [1, 1, 1] => 0.1,
  [1, 1, 2] => 0.2,
  [1, 2, 1] => 0.3,
  [1, 2, 2] => 0.4,
  # Locus 2
  [2, 1, 1] => 0.4,
  [2, 1, 2] => 0.3,
  [2, 2, 1] => 0.2,
  [2, 2, 2] => 0.1,
  # Locus 3
  [3, 1, 1] => 0.2,
  [3, 1, 2] => 0.1,
  [3, 2, 1] => 0.4,
  [3, 2, 2] => 0.3
)

ls = NK.NKLandscape(n, k, a, links, contribs)

pop = [
  1 1 1;
  1 1 2;
  1 2 1;
  1 2 2;
  2 1 1;
  2 1 2;
  2 2 1;
  2 2 2;
] |> transpose

fits = [
  (0.1 + 0.4 + 0.2) / 3.0,
  (0.1 + 0.3 + 0.4) / 3.0,
  (0.2 + 0.2 + 0.2) / 3.0,
  (0.2 + 0.1 + 0.4) / 3.0,
  (0.3 + 0.4 + 0.1) / 3.0,
  (0.3 + 0.3 + 0.3) / 3.0,
  (0.4 + 0.2 + 0.1) / 3.0,
  (0.4 + 0.1 + 0.3) / 3.0,
]

facts("NKLandscapes.jl fast functional tests") do

  context("NK.fitness(...)") do
    function testfitness(g, f0)
      f1 = NK.fitness(g, ls)
      @fact f0 --> roughly(f1) "Expected $(f0), found $(f1), for $(g)"
    end

    for i = 1:8
      testfitness(pop[:,i], fits[i])
    end
  end

  context("NK.popfits(...)") do
    fitsa = NK.popfits(pop, ls)

    for i = 1:8
      @fact fits[i] --> roughly(fitsa[i]) "Expected $(fits[i]), found $(fitsa[i]) for $(i)th"
    end
  end

  context("NK.propsel(...)") do
    srand(0)
    newpop = NK.propsel(pop, ls)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop[:,i] == pop[:,j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop, ls))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.tournsel(...)") do
    srand(0)
    newpop = NK.tournsel(pop, ls, 8)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop[:,i] == pop[:,j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop, ls))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.moransel(...)") do
    srand(0)
    newpop = NK.moransel(pop, ls, 8)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop[:,i] == pop[:,j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop, ls))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.bwmutate(...)") do
    srand(0)
    trials = 1000
    mutprob = 0.1
    total = 0.0
    for _ = 1:trials
      newpop = NK.bwmutate(pop, ls, mutprob)
      total += (pop .!= newpop) |> sum
    end
    meanmuts = total / (trials * length(pop))
    @fact mutprob --> roughly(meanmuts; atol=mutprob / 10) "Expected $(mutprob) mutrate, found $(meanmuts)"
  end

  context("NK.bsmutate(...)") do
    srand(0)
    newpop = NK.bsmutate(pop, ls, 1.0)

    for i = 1:8
      count = 0
      for j = 1:ls.n
        if newpop[j,i] != pop[j,i]
          count += 1
        end
      end
      @fact count --> 1 "Expected 1 difference in $(pop[:,i]), $(newpop[:,i]), found $(count)"
    end
  end
end

