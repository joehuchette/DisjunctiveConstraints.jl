module DisjunctiveConstraints

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

include("sets.jl")
include("activity.jl")
include("big_m.jl")

end # module
