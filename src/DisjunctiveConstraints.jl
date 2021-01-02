module DisjunctiveConstraints

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

export DisjunctiveSet, CombinatorialDisjunctiveSet
export minimum_activity, maximum_activity

include("sets.jl")
include("activity.jl")
include("big_m.jl")

end # module
