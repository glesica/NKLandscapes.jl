pcount = 4
psize = 2

srand(1)

# TODO: Links gets set up properly on instantiation.

context("metapopulation.jl") do
  landscapes = [
    NK.NKLandscape(4, 1),
    NK.NKqLandscape(4, 1, 2),
    NK.NKpLandscape(4, 1, 0.90),
  ]

  for ls = landscapes
    mpr = rand(NK.MetaPopulation, ls, psize, pcount)
    mpz = zeros(NK.MetaPopulation, ls, psize, pcount)

    context("$(ls)") do
      context("popsize(...)") do
        @fact NK.popsize(mpr) --> psize
        @fact NK.popsize(mpz) --> psize
      end

      context("popct(...)") do
        @fact NK.popct(mpr) --> pcount
        @fact NK.popct(mpz) --> pcount
      end

      context("popfits(...)") do
        @fact size(NK.popfits(mpr))[1] --> psize
        @fact size(NK.popfits(mpz))[1] --> psize
        @fact size(NK.popfits(mpr))[2] --> pcount
        @fact size(NK.popfits(mpz))[2] --> pcount
      end
    end
  end
end
