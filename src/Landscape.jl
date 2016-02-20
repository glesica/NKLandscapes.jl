import Base.show

export Landscape, NKLandscape, NKqLandscape, NKpLandscape, prob_neutral

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
  links = zeros(Int64, k + 1, n)
  if near
    if isodd(k)
      error("k must be even to use nearest-neighbor interactions")
    end
    ks = div(k, 2)
    for i = 1:n
      lr = links[:,i]
      lr = [j for j = (i - ks):(i + ks)]
      # Periodic boundary.
      lr[lr .< 1] = lr[lr .< 1] + n
      lr[lr .> n] = lr[lr .> n] - n
      links[:,i] = lr
    end
  else
    for i = 1:n
      links[:,i] = [i; sample(setdiff(1:n, [i]), k, replace=false) |> sort]
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

# TODO: The functions below should move to Enumeration.jl.

@doc """ fitrange(ls::Landscape)

Returns a tuple of the min and max fitness for a landscape. Note that
this only works if `a = 2` and it won't work for large values of `N`.
"""
function fitrange(ls::Landscape)
  fmin = 1.0
  fmax = 0.0
  for g = [[bit == '1' ? 2 : 1 for bit in bin(n, ls.n)] for n in 0:2^ls.n-1]
    f = fitness(g, ls)
    if f < fmin
      fmin = f
    end
    if f > fmax
      fmax = f
    end
  end
  return Float64[fmin, fmax]
end

function leastfit(ls::Landscape)
end

function mostfit(ls::Landscape)
end
