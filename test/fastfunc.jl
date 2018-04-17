srand(1)

n = 3
k = 1
a = 2
# [Locus 1, Locus 2, Locus 3]
links = [0b101, 0b011, 0b110]
contribs = NK.AlleleContribs([
  Dict{NK.AlleleString,Float64}(
    # Locus 1
    0b000 => 0.1,
    0b100 => 0.2,
    0b001 => 0.3,
    0b101 => 0.4
  ),
  Dict{NK.AlleleString,Float64}(
    # Locus 2
    0b000 => 0.4,
    0b001 => 0.3,
    0b010 => 0.2,
    0b011 => 0.1,
  ),
  Dict{NK.AlleleString,Float64}(
    # Locus 3
    0b000 => 0.2,
    0b010 => 0.1,
    0b100 => 0.4,
    0b110 => 0.3
  )
])

ls = NK.NKLandscape(n, k, a, links, contribs)
pop = NK.Population([NK.Genotype(alleles, ls) for alleles = 0:7])

fits = [
  (0.1 + 0.4 + 0.2) / 3.0,
  (0.3 + 0.3 + 0.2) / 3.0,
  (0.1 + 0.2 + 0.1) / 3.0,
  (0.3 + 0.1 + 0.1) / 3.0,
  (0.2 + 0.4 + 0.4) / 3.0,
  (0.4 + 0.3 + 0.4) / 3.0,
  (0.2 + 0.2 + 0.3) / 3.0,
  (0.4 + 0.1 + 0.3) / 3.0
]

context("Fast functional tests") do
  context("NK.fitness(...)") do
    function testfitness(g, f0)
      f1 = NK.fitness(g)
      @fact f0 --> roughly(f1) "Expected $(f0), found $(f1), for $(g)"
    end

    for i = 1:8
      testfitness(pop.genotypes[i], fits[i])
    end
  end

  context("NK.popfits(...)") do
    fitsa = NK.popfits(pop)

    for i = 1:8
      @fact fits[i] --> roughly(fitsa[i]) "Expected $(fits[i]), found $(fitsa[i]) for $(i)th"
    end
  end

  context("NK.propsel(...)") do
    srand(0)
    newpop = NK.propsel(pop)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop.genotypes[i] == pop.genotypes[j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.tournsel(...)") do
    srand(0)
    newpop = NK.tournsel(pop, 8)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop.genotypes[i] == pop.genotypes[j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.moransel(...)") do
    srand(0)
    newpop = NK.moransel(pop, 16)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop.genotypes[i] == pop.genotypes[j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
  end

  context("NK.elitesel(...)") do
    srand(0)
    newpop = NK.elitesel(pop, 2)

    counts = zeros(Int64, 8)
    for i = 1:8
      for j = 1:8
        if newpop.genotypes[i] == pop.genotypes[j]
          counts[j] += 1
        end
      end
    end
    @fact any(c -> c > 1, counts) --> true "Expected a different population after selection"
    @fact newpop.genotypes |> unique |> length --> 2 "Expected two genotypes in new population"
    oldmean = mean(fits)
    newmean = mean(NK.popfits(newpop))
    @fact oldmean --> less_than(newmean) "Expected greater mean fitness after selection"
    elites = sort(pop.genotypes, by=(g) -> NK.fitness(g))[end-1:end]
    for g = newpop.genotypes
      @fact g --> anyof(elites...) "Expected members of new population to be elites"
    end
    for g = elites
      @fact g --> anyof(newpop.genotypes...) "Expected all elites to be present in new population"
    end
  end

  context("NK.bwmutate(...)") do
    srand(0)
    trials = 1000
    mutprob = 0.1
    total = 0.0
    for _ = 1:trials
      newpop = NK.bwmutate(pop, mutprob)
      for i = 1:NK.popsize(pop)
        total += (xor(pop.genotypes[i].alleles, newpop.genotypes[i].alleles)) |> count_ones
      end
    end
    meanmuts = total / (trials * NK.popsize(pop) * n)
    @fact mutprob --> roughly(meanmuts; atol=mutprob / 10) "Expected $(mutprob) mutrate, found $(meanmuts)"
  end

  context("NK.bsmutate(...)") do
    srand(0)
    newpop = NK.bsmutate(pop, 1.0)

    for i = 1:NK.popsize(pop)
      delta = abs(count_ones(newpop.genotypes[i].alleles) - count_ones(pop.genotypes[i].alleles))
      @fact delta --> 1 "Expected 1 difference in $(pop.genotypes[i]), $(newpop.genotypes[i]), found $(delta)"
    end
  end
end

