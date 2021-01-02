abstract type AbstractActivityMethod end

struct IntervalArithmetic <: AbstractActivityMethod end

struct LinearProgrammingRelaxation <: AbstractActivityMethod
    opt_factory::Any
end

struct FullFormulation <: AbstractActivityMethod
    opt_factory::Any
end

struct PartiallyEnumerateDisjunction{T<:Real} <: AbstractActivityMethod
    disjunction_ci::DisjunctionCI{T}
    fixings::Vector{Int}
    sub_method::AbstractActivityMethod
end

function minimum_activity(
    model::MOI.ModelLike,
    x,
    method::AbstractActivityMethod,
)
    return -maximum_activity(model, -x, method)
end

function maximum_activity(
    model::MOI.ModelLike,
    x::T,
    method::AbstractActivityMethod,
) where {T<:Real}
    return maximum_activity(
        model,
        convert(MOI.ScalarAffineFunction{T}, x),
        method,
    )
end

function maximum_activity(
    model::MOI.ModelLike,
    aff::MOI.ScalarAffineFunction{T},
    method::IntervalArithmetic,
)::T where {T<:Real}
    val = aff.constant
    for term in aff.terms
        lb, ub = MOIU.get_bounds(model, T, term.variable_index)
        coeff = term.coefficient
        if coeff > 0
            val += coeff * ub
        elseif coeff < 0
            val += coeff * lb
        end
    end
    return val
end

function _is_linear_constraint(::MOI.ConstraintIndex{F,S}) where {F,S}
    if !(
        F <: Union{
            MOI.SingleVariable,
            MOI.ScalarAffineFunction,
            MOI.VectorOfVariables,
            MOI.VectorAffineFunction,
        }
    )
        return false
    end
    if !(
        S <: Union{
            MOI.EqualTo,
            MOI.LessThan,
            MOI.GreaterThan,
            MOI.Interval,
            MOI.Zeros,
            MOI.Nonnegatives,
            MOI.Nonpositives,
        }
    )
        return false
    end
    return true
end

_get_constraint_filter(::LinearProgrammingRelaxation) = _is_linear_constraint
_get_constraint_filter(::FullFormulation) = nothing

function maximum_activity(
    model::MOI.ModelLike,
    aff::MOI.ScalarAffineFunction{T},
    method::Union{LinearProgrammingRelaxation,FullFormulation},
)::T where {T<:Real}
    # TODO: This extra copy is to handle optimizers that implement their own
    #       MOI.copy_to, and so crash if you try to call MOI.automatic_copy_to.
    #       We really should implement our own copy_to variant.
    placeholder_opt = MOIU.Model{T}()
    MOIU.automatic_copy_to(
        placeholder_opt,
        model,
        copy_names = false,
        filter_constraints = _get_constraint_filter(method),
    )
    MOI.set(
        placeholder_opt,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
        aff,
    )
    MOI.set(placeholder_opt, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    opt = method.opt_factory()
    # TODO: What happens to ZeroOne constraints? Should get relaxed to [0,1]?
    MOI.copy_to(opt, placeholder_opt, copy_names = false)
    MOI.optimize!(opt)
    term = MOI.get(opt, MOI.TerminationStatus())
    if term == MOI.INFEASIBLE
        return T(-Inf)
    elseif term == MOI.DUAL_INFEASIBLE
        return T(Inf)
    elseif term == MOI.OPTIMAL
        return T(MOI.get(opt, MOI.ObjectiveValue()))
    else
        @info "Unusual solution status $term; returning valid dual bound."
        return T(MOI.get(opt, MOI.ObjectiveBound()))
    end
end

function maximum_activity(
    model::MOI.ModelLike,
    aff::MOI.ScalarAffineFunction{T},
    method::PartiallyEnumerateDisjunction,
)::T where {T<:Real}
    vector_f = MOI.get(model, MOI.ConstraintFunction(), method.disjunction_ci)
    vector_s = MOI.get(model, MOI.ConstraintSet(), method.disjunction_ci)
    max_val = T(-Inf)
    for i in method.fixings
        # TODO: Copy once and change constraint right-hand sides. However, this
        # only works for non-ranged constraints whose sense is the same for
        # each alternative.
        constrained_model = MOIU.Model{T}()
        MOI.copy_to(constrained_model, model, copy_names = false)
        for (j, scalar_f) in enumerate(MOIU.scalarize(vector_f))
            lb = vector_s.lbs[i][j]
            ub = vector_s.ubs[i][j]
            scalar_s = MOI.Interval(lb, ub)
            MOI.add_constraint(constrained_model, scalar_f, scalar_s)
        end
        max_val = max(
            max_val,
            maximum_activity(constrained_model, aff, method.sub_method),
        )
    end
    return max_val
end
