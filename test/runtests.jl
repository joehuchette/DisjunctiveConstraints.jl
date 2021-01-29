using Test
using DisjunctiveConstraints
import Cbc, MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

include("sets.jl")
include("activity.jl")
include("big_m.jl")
