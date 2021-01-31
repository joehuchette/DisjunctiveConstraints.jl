abstract type AbstractDisjunctiveFormulation end

struct NaiveBigM <: AbstractDisjunctiveFormulation
    activity_method::AbstractActivityMethod
end

struct NaiveBigMState
    z_vis::Vector{VI}
    sum_ci::CI{SAF,ET}
    lt_cis::Matrix{Union{Nothing,CI{SAF,LT}}}
    gt_cis::Matrix{Union{Nothing,CI{SAF,GT}}}
end

function formulate!(
    model::MOI.ModelLike,
    method::NaiveBigM,
    disjunction::Disjunction,
)
    z_v, z_c = MOI.add_constrained_variables(
        model,
        [ZO() for i in 1:num_alternatives(disjunction)],
    )
    return formulate!(model, method, disjunction, z_v)
end

function formulate!(
    model::MOI.ModelLike,
    method::NaiveBigM,
    disjunction::Disjunction,
    z_vis::Vector{VI},
)
    # We need to copy the model (do we?) because we are intermixing adding
    # constraints and optimizing over the model.
    # TODO: Cache this activity_model inside the method for resuse.
    activity_model = MOIU.Model{Float64}()
    MOI.copy_to(activity_model, model, copy_names = false)
    m = MOI.output_dimension(disjunction.f)
    @assert m ==
            MOI.dimension(disjunction.s) ==
            size(disjunction.s.lbs, 1) ==
            size(disjunction.s.ubs, 1)
    d = length(z_vis)
    @assert d ==
            num_alternatives(disjunction) ==
            size(disjunction.s.lbs, 2) ==
            size(disjunction.s.ubs, 2)

    sum_ci = MOI.add_constraint(
        model,
        SAF([SAT(1.0, vi) for vi in z_vis], 0.0),
        MOI.EqualTo(1.0),
    )

    f_scalar = MOIU.scalarize(disjunction.f)

    lt_cis = similar(disjunction.s.lbs, Union{Nothing,CI{SAF,LT}})
    gt_cis = similar(disjunction.s.ubs, Union{Nothing,CI{SAF,GT}})
    for i in 1:m
        f = f_scalar[i]
        for j in 1:d
            let lb = disjunction.s.lbs[i, j]
                gt_cis[i, j] = (
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
                        m_val = minimum_activity(
                            method.activity_method,
                            activity_model,
                            f,
                        )
                        if m_val == -Inf
                            throw(
                                ValueError(
                                    "Encountered infinite big-M coefficient, cannot formulate.",
                                ),
                            )
                        end
                        MOI.add_constraint(
                            model,
                            f + (lb - m_val) * SV(z_vis[j]),
                            GT(lb),
                        )
                    end
                )
            end

            let ub = disjunction.s.ubs[i, j]
                lt_cis[i, j] = (
                    if ub == -Inf
                        throw(
                            ArgumentError(
                                "The disjunction is trivially infeasible.",
                            ),
                        )
                    elseif ub == Inf
                        # Do nothing
                        nothing
                    else
                        m_val = maximum_activity(
                            method.activity_method,
                            activity_model,
                            f,
                        )
                        if m_val == Inf
                            throw(
                                ValueError(
                                    "Encountered infinite big-M coefficient, cannot formulate.",
                                ),
                            )
                        end
                        MOI.add_constraint(
                            model,
                            f + (ub - m_val) * SV(z_vis[j]),
                            LT(ub),
                        )
                    end
                )
            end
        end
    end
    return NaiveBigMState(z_vis, sum_ci, lt_cis, gt_cis)
end
