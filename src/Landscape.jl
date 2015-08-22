using Distributions

export Landscape, prob_neutral

typealias Contribs Dict{Vector{Int64}, Float64}

type Landscape
  n::Int64
  k::Int64
  a::Int64
  p::Float64
  links::Matrix{Int64}
  contribs::Contribs
end

function Landscape(n::Int64, k::Int64, p::Float64, a::Int64, nearest::Bool=false)
  if k >= n
    error("k must be strictly less than n")
  end
  links = zeros(Int64, k + 1, n)
  if nearest
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
      links[:,i] = [i, sample(setdiff(1:n, [i]), k, replace=false) |> sort]
    end
  end
  contribs = Contribs()
  return Landscape(n, k, a, p, links, contribs)
end

Landscape(n::Int64, k::Int64) = Landscape(n, k, 0.0, 2)
Landscape(n::Int64, k::Int64, p::Float64) = Landscape(n, k, p, 2)
Landscape(n::Int64, k::Int64, a::Int64) = Landscape(n, k, 0.0, a)
Landscape(n::Int64, k::Int64, near::Bool) = Landscape(n, k, 0.0, 2, near)
Landscape(n::Int64, k::Int64, p::Float64, near::Bool) = Landscape(n, k, p, 2, near)
Landscape(n::Int64, k::Int64, a::Int64, near::Bool) = Landscape(n, k, 0.0, a, near)

# Barnett, L. (1998). Ruggedness and Neutrality - The NKp Family of Fitness
# Landscapes. In Artificial Life VI (pp. 18–27).
function prob_neutral(l::Landscape)
  l.p ^ 2 * (1 - l.k / (l.n - 1) * (1 - l.p ^ 2)) ^ (l.n - 1)
end
