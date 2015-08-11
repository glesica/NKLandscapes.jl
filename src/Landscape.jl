using Distributions

export Landscape

type Landscape
  n::Int64
  k::Int64
  a::Int64
  links::Matrix{Int64}
  contribs::Array{Float64}
end

function Landscape(n::Int64, k::Int64, p::Float64, a::Int64)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(Int64, n, k + 1)
  for i = 1:n
    links[i,:] = [i, sample(setdiff(1:n, [i]), k, replace=false) |> sort]
  end
  contribs = rand(Float64, n, (ones(Int, k + 1) * a)...)
  # Zero contributions with probability p.
  contribs[rand(Float64, size(contribs)...) .< p] = 0.0
  return Landscape(n, k, a, links, contribs)
end

Landscape(n::Int64, k::Int64) = Landscape(n, k, 0.0, 2)
Landscape(n::Int64, k::Int64, p::Float64) = Landscape(n, k, p, 2)
Landscape(n::Int64, k::Int64, a::Int64) = Landscape(n, k, 0.0, a)

