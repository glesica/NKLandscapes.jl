export AlleleMask, AlleleString, AlleleLinks, AlleleContribs

@doc """AlleleMask UInt128

A bit mask that identifies one or several alleles within a genotype.
"""
typealias AlleleMask UInt128

@doc """AlleleString UInt128

A bit string that represents the alleles of a particular genotype.
"""
typealias AlleleString UInt128

@doc """AlleleLinks Vector{AlleleMask}

A vector of masks specifying the epistatic links for each allele. The ith
element of the vector stores a mask of alleles that are epistatically linked
to the allele in the ith least significant position.

Each mask should have a 1 in the position that corresponds to the locus
discussed above, and in the K positions that correspond to the loci to which it
is epistatically linked.
"""
typealias AlleleLinks Vector{AlleleMask}

@doc """AlleleContribs Vector{Dict{AlleleString, Float64}}

A data structure in which fitness contributions may be stored. The ith element
of the vector stores contributions of the allele in the ith least significant
position. For example, the 4th element of the vector corresponds to the locus
which contains a 1 bit in the bit string 1000.

Each dictionary then maps a bit string into a fitness contribution. The bit string should contain all zeros except for the following:

  * The locus that corresponds to the dictionary should contain its actual
    value, either 1 or 0
  * The K loci to which the locus above is linked should contain their actual
    values, either 1 or 0
"""
typealias AlleleContribs Vector{Dict{AlleleString, Float64}}

newAlleleContribs(n::Int64) = [Dict{AlleleString, Float64}() for _ = 1:n]

