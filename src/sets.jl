abstract type AbstractDisjunctiveSet <: MOI.AbstractVectorSet end

# For alternative $j$, impose that entry $i$ lies in [lbs[i,j], ubs[i,j]]
struct DisjunctiveSet{T<:Real} <: AbstractDisjunctiveSet
    lbs::Matrix{T}
    ubs::Matrix{T}

    function DisjunctiveSet(
        lbs::Matrix{R},
        ubs::Matrix{S},
    ) where {R<:Real,S<:Real}
        if size(lbs) != size(ubs)
            throw(
                DimensionMismatch(
                    "Dimensions of lowerbounds and upperbounds must match.",
                ),
            )
        end
        return new{promote_type(R, S)}(lbs, ubs)
    end
end

MOI.dimension(ds::DisjunctiveSet) = size(ds.lbs, 1)
num_alternatives(ds::DisjunctiveSet) = size(ds.lbs, 2)

const DisjunctionCI{T} =
    MOI.ConstraintIndex{MOI.VectorAffineFunction{T},DisjunctiveSet}

# Constrains f in Simplex AND sparsity matches disjunction
struct CombinatorialDisjunctiveSet <: AbstractDisjunctiveSet
    sparsity::Vector{Vector{Int}}
    dim::Int

    function CombinatorialDisjunctiveSet(sparsity::Vector{Vector{Int}})
        dim = 0
        for disjunct in sparsity
            if minimum(disjunct) <= 0
                throw(ArgumentError("Nonpositive variable index."))
            end
            dim = max(0, maximum(disjunct))
        end
        return new(sparsity, dim)
    end
end

MOI.dimension(cdc::CombinatorialDisjunctiveSet) = cdc.dim

struct Disjunction{F<:MOI.AbstractVectorFunction,S<:AbstractDisjunctiveSet}
    f::F
    s::S

    function Disjunction{F,S}(
        f::F,
        s::S,
    ) where {F<:MOI.AbstractVectorFunction,S<:AbstractDisjunctiveSet}
        @assert dimension(f) == dimension(s)
        return new(f, s)
    end
end

num_alternatives(disj::Disjunction) = num_alternatives(disj.s)
