using NKLandscapes
using Base.Test

const BOUNDARIES = ["Nearest-neighbor", "Random"]
const LC = 100 # Landscape count
const PAIRS = [(8, 0), (8, 2), (8, 4), (8, 8),
(16, 0), (16, 2), (16, 4), (16, 8), (16, 16),
(24, 0), (24, 2), (24, 4), (24, 8), (24, 16), (24, 24),
(48, 0), (48, 2), (48, 4), (48, 8), (48, 16), (48, 24), (48, 48),
(96, 0), (96, 2), (96, 4), (96, 8), (96, 16), (96, 24), (96, 48), (96, 96)]

# Mean fitness of local optima (nearest-neighbor interactions) & Mean walk
# lengths to local optima (nearest-neighbor interactions) Kauffman, p. 55-56

# Mean fitness of local optima (random interactions) & Mean walk lengths to
# local optima (random interactions) Kauffman, p. 57

println("Length of Paths to Optima and Fitness of Optima")
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
      l = NKLandscape(n, k, near=b == 1)
      g = rand(Genotype, l)
      w = random_adaptive_walk(g, l)
      lengths[i] = w.length
      fitnesses[i] = w.fitnesses[end]
    end
    @printf("%d\t%d\t%.2f\t%.2f\t%.3f\t%.3f\n", n, k, mean(lengths), std(lengths), mean(fitnesses), std(fitnesses))
  end
end

# Number of optima. Kauffman, p. 60.

MAX_FOUND = 10_000
MAX_FAILS = 100

println("Number of Optima")
for b = [1, 2]
  println(BOUNDARIES[b], " Interactions")
  println("N\tK\tμ\tσ")
  for (n, k) = PAIRS[1:9]  # Only runs tests for N = 8 and N = 16.  Modify to run more tests (time consuming).
    if n == k
      k = k - 1
      if b == 1
        k = k - 1
      end
    end
    counts = zeros(LC)
    for i = 1:LC
      since_last = 0
      found_opts = Set()
      l = NKLandscape(n, k, near=b == 1)
      while length(found_opts) <= MAX_FOUND && since_last <= MAX_FAILS
        g = rand(Genotype, l)
        w = random_adaptive_walk(g, l)
        opt = w.history_list[:,end]
        if !(opt in found_opts)
          push!(found_opts, opt)
          since_last = 0
        else
          since_last += 1
        end
      end
      counts[i] = length(found_opts)
    end
    @printf("%d\t%d\t%.2f\t%.2f\n", n, k, mean(counts), std(counts))
  end
end


