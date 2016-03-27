import Base.show

export Landscape, prob_neutral

@doc """ function makelinks(n::Int64, k::Int64, near::Bool)

Return an array of links (bit masks) that select the linked loci.
Thus, link[i] selects the loci that are linked to locus  `i`.
Each link contains `k + 1` one bits.  link[i] is guaranteed to have bit i 
(counting from the right) set to one.  The positions of the remaining `k` 
one bits are randomly selected from the set {1:n} \ {i} (where \ denotes
set difference).
"""
function makelinks(n::Int64, k::Int64, near::Bool)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(AlleleMask, n)
  if near
    # TODO: See if we can come up with a prettier algorithm here.
    if isodd(k)
      error("Periodic boundary requires an even-valued K.")
    end
    offset = div(k, 2)
    for i = 1:n
      mask = AlleleMask(0)
      for j = (i - offset):(i + offset)
        shift = if j < 1
          j + n - 1
        elseif j > n
          j - n - 1
        else
          j - 1
        end
        mask = mask | (AlleleMask(1) << shift)
      end
      links[i] = mask
    end
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

