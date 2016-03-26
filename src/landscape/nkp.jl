export NKpLandscape

# Construct an NKp landscape, see:
# Barnett, L. (1998). Ruggedness and Neutrality - The NKp Family of Fitness
# Landscapes. In Artificial Life VI (pp. 18–27).
type NKpLandscape <: Landscape
  n::Int64
  k::Int64
  p::Float64
  a::Int64
  links::AlleleLinks
  contribs::AlleleContribs
end

function NKpLandscape(n::Int64, k::Int64, p::Float64; a::Int64=2, near::Bool=false)
  if a != 2
    # TODO:  implement a values other than 2
    error("only a == 2 is implemented at this time")
  end
  links = makelinks(n, k, near)
  contribs = newAlleleContribs(n)
  return NKpLandscape(n, k, p, a, links, contribs)
end

show(io::Base.IO, l::NKpLandscape) = print(io, "NKpLandscape($(l.n), $(l.k), $(l.p))")

# Compute an estimate of the probability that a mutation within an NKp
# landscape is neutral.
# Barnett, L. (1998). Ruggedness and Neutrality - The NKp Family of Fitness
# Landscapes. In Artificial Life VI (pp. 18–27).
function prob_neutral(l::NKpLandscape)
  l.p ^ 2 * (1 - l.k / (l.n - 1) * (1 - l.p ^ 2)) ^ (l.n - 1)
end

