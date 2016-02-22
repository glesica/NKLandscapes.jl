export NKLandscape

# Construct an NK lanscape, see:
# Kauffman, S. A. (1993). The Origins of Order: Self-Organization and Selection
# in Evolution. New York, New York, USA: Oxford University Press.
type NKLandscape <: Landscape
  n::Int64
  k::Int64
  a::Int64
  links::AlleleLinks
  contribs::AlleleContribs
end

function NKLandscape(n::Int64, k::Int64; a::Int64=2, near::Bool=false)
  links = makelinks(n, k, near)
  contribs = AlleleContribs()
  return NKLandscape(n, k, a, links, contribs)
end

show(io::Base.IO, l::NKLandscape) = print(io, "NKLandscape($(l.n), $(l.k))")

