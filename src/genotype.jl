import Base.Random: rand, zeros
import Base: ==, hash, show

export Genotype, contribs, fitness

@doc """A genotype representation.
"""
type Genotype{T <: Landscape}
  alleles::AlleleString
  landscape::T
end

Genotype{T <: Landscape}(alleles::Integer, landscape::T) = Genotype(AlleleString(alleles), landscape)

Genotype{T <: Landscape}(g::Genotype{T}) = Genotype(g.alleles, g.landscape)

@doc """contribs(g::Genotype, update::Function)

Return a vector of contributions where the ith element in
the vector is the contribution made by the ith allele in the
genotype, given the values of the k alleles to which it is
epistatically linked.
"""
function contribs(g::Genotype, update::Function)
  return map(1:g.landscape.n) do i
    linksmask = AlleleMask(g.landscape.links[i])
    contribstring = AlleleString(g.alleles & linksmask)
    return get!(update, g.landscape.contribs[i], contribstring)
  end
end

@doc """contribs(g::Genotype)
"""
contribs(g::Genotype{NKLandscape}) = contribs(g, () -> rand())

@doc """contribs(g::Genotype)
"""
contribs(g::Genotype{NKqLandscape}) = contribs(g, () -> rand(0:(g.landscape.q - 1)))

@doc """contribs(g::Genotype)
"""
function contribs(g::Genotype{NKpLandscape})
  update = () -> begin
    if rand() < g.landscape.p
      return 0.0
    else
      return rand()
    end
  end
  return contribs(g, update)
end

# TODO: Do we include the fake NKp zeros in the sum or just let the fitness range lower?
@doc """fitness(g::Genotype{NKqLandscape})
"""
fitness(g::Genotype{NKqLandscape}) = mean(contribs(g)) / (g.landscape.q - 1)

@doc """fitness(g::Genotype)

Compute the fitness of a particular genotype.
"""
fitness(g::Genotype) = mean(contribs(g))

function rand(::Type{Genotype}, ls::Landscape)
  nmask = (AlleleMask(1) << ls.n) - 1
  Genotype(rand(AlleleString) & nmask, ls)
end

zeros(::Type{Genotype}, ls::Landscape) = Genotype(0, ls)

==(left::Genotype, right::Genotype) = left.alleles == right.alleles && ===(left.landscape, right.landscape)

hash(g::Genotype) = xor(hash(Genotype), xor(hash(g.landscape), hash(g.alleles)))

show(io::Base.IO, g::Genotype) = print(io, "Genotype($(bits(g.alleles)[(end - g.landscape.n + 1):end]), $(g.landscape))")

