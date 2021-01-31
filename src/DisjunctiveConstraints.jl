module DisjunctiveConstraints

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex
const SV = MOI.SingleVariable
const SAT = MOI.ScalarAffineTerm{Float64}
const SAF = MOI.ScalarAffineFunction{Float64}
const ET = MOI.EqualTo{Float64}
const LT = MOI.LessThan{Float64}
const GT = MOI.GreaterThan{Float64}
const ZO = MOI.ZeroOne

export DisjunctiveSet, CombinatorialDisjunctiveSet
export minimum_activity, maximum_activity

include("sets.jl")
include("activity.jl")
include("big_m.jl")

end # module
