module NKLandscapes

include("Landscape.jl")
include("Genotype.jl")
include("Population.jl")
include("MetaPopulation.jl")
include("Neighbors.jl")
include("Walks.jl")
include("NeutralNets.jl")

include("selection/moran.jl")
include("selection/proportional.jl")
include("selection/tournament.jl")

include("Mutation.jl")

include("Migration.jl")

end

