import NKLandscapes
const NK = NKLandscapes

n = 3
k = 1
a = 2
links = [1 2 3;
         2 3 1]
contribs = NK.Contribs(
  # Locus 1
  [1, 1, 1] => 0.1,
  [1, 1, 2] => 0.2,
  [1, 2, 1] => 0.3,
  [1, 2, 2] => 0.4,
  # Locus 2
  [2, 1, 1] => 0.4,
  [2, 1, 2] => 0.3,
  [2, 2, 1] => 0.2,
  [2, 2, 2] => 0.1,
  # Locus 3
  [3, 1, 1] => 0.2,
  [3, 1, 2] => 0.1,
  [3, 2, 1] => 0.4,
  [3, 2, 2] => 0.3
)

ls = NK.NKLandscape(n, k, a, links, contribs)

pop = [
  1 1 1;
  1 1 2;
  1 2 1;
  1 2 2;
  2 1 1;
  2 1 2;
  2 2 1;
  2 2 2;
] |> transpose

fits = [
  (0.1 + 0.4 + 0.2) / 3.0,
  (0.1 + 0.3 + 0.4) / 3.0,
  (0.2 + 0.2 + 0.2) / 3.0,
  (0.2 + 0.1 + 0.4) / 3.0,
  (0.3 + 0.4 + 0.1) / 3.0,
  (0.3 + 0.3 + 0.3) / 3.0,
  (0.4 + 0.2 + 0.1) / 3.0,
  (0.4 + 0.1 + 0.3) / 3.0,
]

# NK.fitness(...)

function testfitness(g, f0)
  f1 = NK.fitness(g, ls)
  @assert isapprox(f0, f1) "Expected $(f0), found $(f1), for $(g)"
end

for i = 1:8
  testfitness(pop[:,i], fits[i])
end

# NK.popfits(...)

fitsa = NK.popfits(pop, ls)

for i = 1:8
  @assert isapprox(fits[i], fitsa[i]) "Expected $(fits[i]), found $(fitsa[i]) for $(i)th"
end

# NK.propsel(...)

srand(0)
newpop = NK.propsel(pop, ls)

counts = zeros(Int64, 8)
for i = 1:8
  for j = 1:8
    if newpop[:,i] == pop[:,j]
      counts[j] += 1
    end
  end
end
@assert any(c -> c != 1, counts) "Expected a different population after selection"

# NK.bwmutate(...)

srand(0)
trials = 1000
mutprob = 0.1
total = 0.0
for _ = 1:trials
  newpop = NK.bwmutate(pop, ls, mutprob)
  total += (pop .!= newpop) |> sum
end
meanmuts = total / (trials * length(pop))
@assert isapprox(mutprob, meanmuts, atol=mutprob / 10) "Expected $(mutprob) mutation rate, found $(meanmuts)"

# NK.bsmutate(...)

srand(0)
newpop = NK.bsmutate(pop, ls)

for i = 1:8
  count = 0
  for j = 1:ls.n
    if newpop[j,i] != pop[j,i]
      count += 1
    end
  end
  @assert count == 1 "Expected 1 difference in $(pop[:,i]), $(newpop[:,i]), found $(count)"
end

