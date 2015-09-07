using Distributions

export Landscape, NKLandscape, NKqLandscape, NKpLandscape, prob_neutral

typealias Links Matrix{Int64}

# Generate a table describing the epistatic links between loci. If `near` is
# `true`, then links will be generated between a locus `i` and its `k / 2`
# neighbors to either side, assuming a periodic boundary at the ends of the
# array.
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

typealias Contribs Dict{Vector{Int64}, Float64}

# TODO: Provide a function to pre-fill contribs (for use with small k).

abstract Landscape

type NKLandscape <: Landscape
  n::Int64
  k::Int64
  a::Int64
  links::Links
  contribs::Contribs
end

# Construct an NK lanscape, see:
# Kauffman, S. A. (1993). The Origins of Order: Self-Organization and Selection
# in Evolution. New York, New York, USA: Oxford University Press.
function NKLandscape(n::Int64, k::Int64; a::Int64=2, near::Bool=false)
  links = makelinks(n, k, near)
  contribs = Contribs()
  return NKLandscape(n, k, a, links, contribs)
end

type NKqLandscape <: Landscape
  n::Int64
  k::Int64
  q::Int64
  a::Int64
  links::Links
  contribs::Contribs
end

# Construct an NKq landscape, see:
# Newman, M., & Engelhardt, R. (1998). Effects of selective neutrality on the
# evolution of molecular species. Proceedings of the Royal Society of London.
# Series B: Biological Sciences, 265(1403)
function NKqLandscape(n::Int64, k::Int64, q::Int64; a::Int64=2, near::Bool=false)
  links = makelinks(n, k, near)
  contribs = Contribs()
  return NKqLandscape(n, k, q, a, links, contribs)
end

type NKpLandscape <: Landscape
  n::Int64
  k::Int64
  p::Float64
  a::Int64
  links::Links
  contribs::Contribs
end

# Construct an NKp landscape, see:
# Barnett, L. (1998). Ruggedness and Neutrality - The NKp Family of Fitness
# Landscapes. In Artificial Life VI (pp. 18–27).
function NKpLandscape(n::Int64, k::Int64, p::Float64; a::Int64=2, near::Bool=false)
  links = makelinks(n, k, near)
  contribs = Contribs()
  return NKpLandscape(n, k, p, a, links, contribs)
end

# Compute an estimate of the probability that a mutation within an NKp
# landscape is neutral.
# Barnett, L. (1998). Ruggedness and Neutrality - The NKp Family of Fitness
# Landscapes. In Artificial Life VI (pp. 18–27).
function prob_neutral(l::NKpLandscape)
  l.p ^ 2 * (1 - l.k / (l.n - 1) * (1 - l.p ^ 2)) ^ (l.n - 1)
end
