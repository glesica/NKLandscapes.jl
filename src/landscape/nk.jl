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
  if a != 2
    # TODO:  implement a values other than 2
    error("only a == 2 is implemented at this time")
  end
  links = makelinks(n, k, near)
  contribs = newAlleleContribs(n)
  return NKLandscape(n, k, a, links, contribs)
end

show(io::Base.IO, l::NKLandscape) = print(io, "NKLandscape($(l.n), $(l.k))")

