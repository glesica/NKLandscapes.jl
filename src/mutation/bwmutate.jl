function bwmutate!(g::Genotype, mutprob::Float64)
  mask = AlleleMask(0)
  current = AlleleMask(1)
  for _ = 1:g.landscape.n
    if rand() < mutprob
      mask = mask | current
    end
    current = current << 1
  end
  g.alleles = g.alleles $ mask
end

