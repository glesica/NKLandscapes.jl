export elitesel, elitesel!

# TODO: Implement for `MetaPopulation`.

@doc """elitesel(p::Population, q::Int64=1)

Create a new population through elite selection, which keeps only the `q`
fittest members of the population and copies them to create the new population.

If `q > 1`, the elite individuals will be copied proportionally based on their
fitnesses. So, as an example, assume the fittest individual has fitness 0.75,
the second fittest individual has fitness 0.5, and the population size is 100.
The fittest individual will be copied 60 times and the other will be copied 40
times to make the new population.

Each elite individual will always show up in the resulting population at least
once, even when its fitness is so low that its share in the population would
round to zero.
"""
function elitesel(p::Population, q::Int64=1)
  np = Population(p)
  elitesel!(np, q)
  return np
end

@doc """elitesel!(p::Population, q::Int64=1)

Conduct elite selection in-place.
"""
function elitesel!(p::Population, q::Int64=1)
  elites = sort(p.genotypes, by=(g) -> fitness(g))[end-q+1:end]
  fits = [fitness(g) for g = elites]
  counts = [round(Int64, f / sum(fits) * popsize(p) |> ceil) for f = fits]
  # We need to adjust counts to achieve the same population size, we do this by
  # penalizing the most-fit individual. This is a naive approach, but it should
  # work fine for most use cases. It would be better to ensure that the
  # most-fit individual always ends up more common than the next most-fit, etc.
  while sum(counts) > popsize(p)
    counts[end] -= 1
  end
  next = 1
  for i = 1:length(elites)
    for _ = 1:counts[i]
      p.genotypes[next] = Genotype(elites[i])
      next += 1
    end
  end
end

