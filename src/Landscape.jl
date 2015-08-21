using Distributions

export Landscape

type Landscape
  n::Int64
  k::Int64
  a::Int64
  p::Float64
  links::Matrix{Int64}
  contribs::Array{Float64}
end

function Landscape(n::Int64, k::Int64, p::Float64, a::Int64, nearest::Bool=false)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(Int64, n, k + 1)
  if nearest
    if isodd(k)
      error("k must be even to use nearest-neighbor interactions")
    end
    ks = div(k, 2)
    for i = 1:n
      lr = links[i,:]
      lr = [j for j = (i - ks):(i + ks)]
      # Periodic boundary.
      lr[lr .< 1] = lr[lr .< 1] + n
      lr[lr .> n] = lr[lr .> n] - n
      links[i,:] = lr
    end
  else
    for i = 1:n
      links[i,:] = [i, sample(setdiff(1:n, [i]), k, replace=false) |> sort]
    end
  end
  contribs = rand(Float64, n, (ones(Int64, k + 1) * a)...)
  # Zero contributions with probability p.
  contribs[rand(Float64, size(contribs)...) .< p] = 0.0
  return Landscape(n, k, a, p, links, contribs)
end

Landscape(n::Int64, k::Int64) = Landscape(n, k, 0.0, 2)
Landscape(n::Int64, k::Int64, p::Float64) = Landscape(n, k, p, 2)
Landscape(n::Int64, k::Int64, a::Int64) = Landscape(n, k, 0.0, a)
Landscape(n::Int64, k::Int64, near::Bool) = Landscape(n, k, 0.0, 2, near)
Landscape(n::Int64, k::Int64, p::Float64, near::Bool) = Landscape(n, k, p, 2, near)
Landscape(n::Int64, k::Int64, a::Int64, near::Bool) = Landscape(n, k, 0.0, a, near)

