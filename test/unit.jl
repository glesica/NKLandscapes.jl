using FactCheck

import NKLandscapes
const NK = NKLandscapes

facts("NKLandscapes") do
  include("genotype.jl")
  include("landscape.jl")
  include("metapopulation.jl")
  include("migration.jl")
  include("mutation.jl")
  include("neighbors.jl")
  include("population.jl")
  include("walks.jl")
  include("basins.jl")

  include("fastfunc.jl")
  include("fastnets.jl")
end

