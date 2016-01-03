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

