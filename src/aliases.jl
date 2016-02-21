export AlleleMask, AlleleString, AlleleLinks, Contribs

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
"""
# FIXME: This only accounts for the case where the allele in question is 1.
typealias AlleleLinks Vector{AlleleMask}

@doc """AlleleContribs Vector{Dict{UInt128, Float64}}

A data structure in which fitness contributions may be stored.  The ith element
of the vector stores contributions of the allele in the ith least significant
position. For example, the 4th element of the vector corresponds to the 1 bit
in the bit string 1000.

Each dictionary then maps a bit string into a fitness contribution.  The bit
string should have 1 bits in the K positions that correspond to the loci
epistatically linked to the locus that corresponds to the vector element.
"""
# FIXME: This only accounts for the case where the allele in question is 1.
typealias AlleleContribs Vector{Dict{AlleleMask, Float64}}

