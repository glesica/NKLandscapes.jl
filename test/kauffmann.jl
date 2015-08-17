using NK
using Base.Test

const LC = 100 # Landscape count
const PAIRS = [(8, 2), (8, 4), (8, 8),
(16, 2), (16, 4), (16, 8), (16, 16)]
#(24, 2), (24, 4), (24, 8), (24, 16), (24, 24),
#(48, 2), (48, 4), (48, 8), (48, 16), (48, 24), (48, 48),
#(96, 2), (96, 4), (96, 8), (96, 16), (96, 24), (96, 48), (96, 96)]

# Mean walk lengths to local optima (random interactions)
# Kauffmann, p. 57

println("N\tK\tμ\tσ")
for (n, k) = PAIRS
  # TODO: Do we actually need to do this?
  if n == k
    k = k - 1
  end
  lengths = zeros(LC)
  for i = 1:LC
    l = Landscape(n, k)
    g0 = random_genome(l)
    f0 = fitness(g0, l)
    while true
      g1 = best_neighbor(g0, l)
      f1 = fitness(g1, l)
      if f1 > f0
        g0 = g1
        f0 = f1
        lengths[i] += 1
      else
        break
      end
    end
  end
  @printf("%d\t%d\t%.2f\t%.2f\n", n, k, mean(lengths), std(lengths))
end

