module DisjunctiveConstraints

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

# For alternative $i$, impose that entry $j$ lies in [lbs[i][j], ubs[i][j]]
struct DisjunctiveSet{T <: Real} <: MOI.AbstractVectorSet
    lbs::Vector{Vector{T}}
    ubs::Vector{Vector{T}}

    function DisjunctiveSet(lbs::Vector{Vector{T}}, ubs::Vector{Vector{T}}) where {T <: Real}
        d = length(lbs)
        if d != length(ubs)
            throw(DimensionMismatch("Inconsistent number of alternatives."))
        end
        if d > 0
            n = length(lbs[1])
            for i in 1:d
                if length(lbs[i]) != n || length(ubs[i]) != n
                    throw(DimensionMismatch("Inconsistent dimension of bounds in alternative $i."))
                end
            end
        end
        return new{T}(lbs, ubs)
    end
end

# Constrains f in Simplex AND sparsity matches disjunction
struct CombinatorialDisjunctiveSet <: MOI.AbstractVectorSet
    sparsity::Vector{Vector{Int}}
    dim::Int

    function CombinatorialDisjunctiveSet(sparsity::Vector{Vector{Int}})
        dim = 0
        for disjunct in sparsity
            if min(disjunct) <= 0
                throw(ArgumentError("Nonpositive variable index."))
            end
            dim = max(0, max(disjunct))
        end
        return new(sparsity, dim)
    end
end

MOI.dimension(ds::DisjunctiveSet) = length(lbs)

const DisjunctionCI{T} = MOI.ConstraintIndex{MOI.VectorAffineFunction{T},DisjunctiveSet}

include("activity.jl")

end # module
