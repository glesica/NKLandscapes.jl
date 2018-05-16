import NKLandscapes
const NK = NKLandscapes

const N = 6
const K = 1
const P = 6  # Population size
const G = 2  # Generations
const E = 2    # Elite carryover
const C = 5    # Intelligence choices
const M = 5  # Moran selection rounds
const T = 1   # Number of trials

out = open("N=$(N)_K=$(K)_P=$(P)_G=$(G)_E=$(E)_C=$(C)_M=$(M)_T=$(T).csv", "w")

function writeheader()
  write(out, join([
    "trial",
    "simulationType",
    "generation",
    "meanFitness",
    "medianFitness",
    "maxFitness"
  ], ","), "\n")
end

function outputrow(trial, simtype, generation, fits)
  write(out, join([
    trial,
    simtype,
    generation,
    mean(fits),
    median(fits),
    maximum(fits)
  ], ","), "\n")
end

writeheader()

for trial = 1:T
  l = NK.NKLandscape(N, K)
  p = rand(NK.Population, l, P)

  # Social
  sp = NK.Population(p)
  for i = 1:G
    fits = NK.popfits(sp)
    outputrow(trial, "social", i, fits)

    NK.bsmutate!(sp, 1.0)
    NK.elitesel!(sp, E)
  end

  # Intelligence
  ip = NK.Population(p)
  for i = 1:G
    fits = NK.popfits(ip)
    outputrow(trial, "intelligence", i, fits)

    for i = 1:NK.popsize(ip)
      g = ip.genotypes[i]
      choice = g
      for _ = 1:C
        option = NK.bsmutate(g, 1.0)
        if NK.fitness(option) > NK.fitness(choice)
          choice = option
        end
      end
      ip.genotypes[i] = choice
    end
    NK.moransel!(ip, M)
  end
end

close(out)

