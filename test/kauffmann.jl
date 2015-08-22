using NK
using Base.Test

const BOUNDARIES = ["Nearest-neighbor", "Random"]
const LC = 100 # Landscape count
const PAIRS = [(8, 0), (8, 2), (8, 4), (8, 8),
(16, 0), (16, 2), (16, 4), (16, 8), (16, 16),
(24, 0), (24, 2), (24, 4), (24, 8), (24, 16), (24, 24),
(48, 0), (48, 2), (48, 4), (48, 8), (48, 16), (48, 24), (48, 48),
(96, 0), (96, 2), (96, 4), (96, 8), (96, 16), (96, 24), (96, 48), (96, 96)]

# Mean fitness of local optima (nearest-neighbor interactions) & Mean walk
# lengths to local optima (nearest-neighbor interactions) Kauffmann, p. 55-56

# Mean fitness of local optima (random interactions) & Mean walk lengths to
# local optima (random interactions) Kauffmann, p. 57

for b = [1, 2]
  println(BOUNDARIES[b], " Interactions")
  println("N\tK\tμ len\tσ len\tμ fit\tσ fit")
  for (n, k) = PAIRS
    if n == k
      k = k - 1
      # Periodic boundary requires even K
      if b == 1
        k = k - 1
      end
    end
    lengths = zeros(LC)
    fitnesses = zeros(LC)
    for i = 1:LC
      l = Landscape(n, k, b == 1)
      g0 = random_genome(l)
      f0 = fitness(g0, l)
      while true
        nbrs = better_neighbors(g0, l)
        if length(nbrs) == 0
          fitnesses[i] = f0
          break
        end
        g1 = nbrs[:,rand(1:end)]
        f1 = fitness(g1, l)
        if f1 > f0
          g0 = g1
          f0 = f1
          lengths[i] += 1
        end
      end
    end
    @printf("%d\t%d\t%.2f\t%.2f\t%.3f\t%.3f\n", n, k, mean(lengths), std(lengths), mean(fitnesses), std(fitnesses))
  end
end

