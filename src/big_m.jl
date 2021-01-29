abstract type AbstractDisjunctiveFormulation end

struct NaiveBigM <: AbstractDisjunctiveFormulation
    activity_method::AbstractActivityMethod
end

function naive_big_m_formulation!(
    model::MOI.ModelLike,
    method::NaiveBigM,
    disjunction::Disjunction,
)
    z_v, z_c = MOI.add_constrained_variables(
        model,
        [MOI.ZeroOne() for i in 1:num_alternatives(disjunction)],
    )
    return naive_big_m_formulation!(model, method, disjunction, z_v)
end

function naive_big_m_formulation!(
    model::MOI.ModelLike,
    method::NaiveBigM,
    disjunction::Disjunction,
    z_vis::Vector{MOI.VariableIndex},
)
    d = dimension(disjunction.f)
    @assert d ==
            dimension(disjunction.s) ==
            length(z_vis) ==
            length(active_mask) ==
            size(disjunction.s.lbs, 2)
    m = size(disjunction.s.lbs, 1)

    sum_ci = MOI.add_constraint(model, MOIU.sum(z_vis), MOI.EqualTo(1.0))

    lt_cis = similar(
        disjunction.s.lbs,
        Union{
            Nothing,
            MOI.ConstraintIndex{
                MOI.ScalarAffineFunction{Float64},
                MOI.LessThan{Float64},
            },
        },
    )
    gt_cis = similar(
        disjunction.s.ubs,
        Union{
            Nothing,
            MOI.ConstraintIndex{
                MOI.ScalarAffineFunction{Float64},
                MOI.GreaterThan{Float64},
            },
        },
    )
    for i in 1:m
        f = disjunction.f[i]
        for j in 1:d
            let lb = disjunction.lbs[j, i]
                gt_cis[j, i] = (
                    if lb == Inf
                        throw(
                            ArgumentError(
                                "The disjunction is trivially infeasible.",
                            ),
                        )
                    elseif lb == -Inf
                        # Do nothing
                        nothing
                    else
                        m_val = minimum_activity(model, f, method.activity_method)
                        if m_val == -Inf
                            throw(
                                ValueError(
                                    "Encountered infinite big-M coefficient, cannot formulate.",
                                ),
                            )
                        end
                        MOI.add_constraint(
                            model,
                            f + (lb - m_val) * z_vis[j],
                            GT(lb),
                        )
                    end
                )
            end

            let ub = disjunction.ubs[j, i]
                lt_cis[j, i] = (
                    if ub == -Inf
                        throw(
                            ArgumentError(
                                "The disjunction is trivially infeasible.",
                            ),
                        )
                    elseif ub == -Inf
                        # Do nothing
                        nothing
                    else
                        m_val = maximum_activity(model, f, method.activity_method)
                        if m_val == Inf
                            throw(
                                ValueError(
                                    "Encountered infinite big-M coefficient, cannot formulate.",
                                ),
                            )
                        end
                        MOI.add_constraint(
                            model,
                            f + (ub - m_val) * z_vis[j],
                            LT(ub),
                        )
                    end
                )
            end
        end
    end
    return sum_ci, lt_cis, gt_cis
end
