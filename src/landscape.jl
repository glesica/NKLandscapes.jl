import Base.show

export Landscape, prob_neutral

# Generate a table describing the epistatic links between loci. If `near` is
# `true`, then links will be generated between a locus `i` and its `k / 2`
# neighbors to either side, assuming a periodic boundary at the ends of the
# array.
#
# Each set of links is a column within the resulting array. There is one
# column for each locus in the genotype. The columns are actually `k + 1`
# long (there are `k + 1` rows in the matrix) because the first row
# points back to the locus whose links are defined by that column. So, the
# first row is just the vector `[1..n]`.
function makelinks(n::Int64, k::Int64, near::Bool)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(AlleleMask, k + 1)
  if near
    # TODO: Periodic boundary
    error("Periodic boundary has been yet been implemented")
  else
    for i = 1:n
      linkedloci = sample(setdiff(0:(n - 1), [i - 1]), k, replace=false)
      mask = reduce((a, b) -> a | (AlleleMask(1) << b),
          AlleleMask(1) << (i - 1), linkedloci)
      links[i] = mask
    end
  end
  return links
end

# TODO: Provide a function to pre-fill contribs (for use with small k).

abstract Landscape

show(io::Base.IO, l::Landscape) = print(io, "$(typeof(l))($(l.n), $(l.k))")

include("landscape/nk.jl")
include("landscape/nkq.jl")
include("landscape/nkp.jl")

