module NKLandscapes

include("Landscape.jl")
include("Genotype.jl")
include("Population.jl")
include("Neighbors.jl")
include("Walks.jl")

include("selection/moran.jl")
include("selection/proportional.jl")
include("selection/tournament.jl")

include("Mutation.jl")

end

