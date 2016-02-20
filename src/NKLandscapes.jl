module NKLandscapes

include("aliases.jl")

include("Landscape.jl")
include("Genotype.jl")
include("Population.jl")
include("MetaPopulation.jl")
include("Neighbors.jl")
include("Walks.jl")
include("Enumeration.jl")

include("selection/moran.jl")
include("selection/proportional.jl")
include("selection/tournament.jl")

include("Mutation.jl")

include("Migration.jl")

end
#using NKLandscapes
