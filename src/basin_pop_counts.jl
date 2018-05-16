using DataStructures

function basin_pop_counts(basins, basinlists, pop::Population)
  pop_counts = copy(basinlists)
  for g in pop.genotypes
    fr = find_root(basins, g.alleles+1)
    at_peak = convert(Int64,g.alleles) == fr-1
    println("fr: ",fr,"  at_peak: ",at_peak)
  end
end

function bwm(g,mp)
    mg = Genotype(copy(g.alleles),g.landscape)
    bwmutate!(mg,mp)
    #println("mg:",mg)
    #println("rg:",rg)
    mg
end

hd( g1, g2) = count_ones( xor(g1.alleles, g2.alleles) )
