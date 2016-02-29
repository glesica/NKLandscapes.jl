export NKqLandscape

# Construct an NKq landscape, see:
# Newman, M., & Engelhardt, R. (1998). Effects of selective neutrality on the
# evolution of molecular species. Proceedings of the Royal Society of London.
# Series B: Biological Sciences, 265(1403)
type NKqLandscape <: Landscape
  n::Int64
  k::Int64
  q::Int64
  a::Int64
  links::AlleleLinks
  contribs::AlleleContribs
end

function NKqLandscape(n::Int64, k::Int64, q::Int64; a::Int64=2, near::Bool=false)
  if a != 2
    # TODO:  implement a values other than 2
    error("only a == 2 is implemented at this time")
  end
  links = makelinks(n, k, near)
  contribs = AlleleContribs(n)
  return NKqLandscape(n, k, q, a, links, contribs)
end

show(io::Base.IO, l::NKqLandscape) = print(io, "NKqLandscape($(l.n), $(l.k), $(l.q))")

