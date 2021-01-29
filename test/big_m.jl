const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex
const SV = MOI.SingleVariable
const SAT = MOI.ScalarAffineTerm{Float64}
const SAF = MOI.ScalarAffineFunction{Float64}
const ET = MOI.EqualTo{Float64}
const LT = MOI.LessThan{Float64}
const GT = MOI.GreaterThan{Float64}
const IN = MOI.Interval{Float64}
const ZO = MOI.ZeroOne

function _build_base_model()
    # Model domain: x in [-1,1]^3 satisfying
    #               x_1 + x_2 <= 0.5 AND
    #               x_1 - x_2 <= 0.6
    model = MOIU.Model{Float64}()
    x_v, x_c =
        MOI.add_constrained_variables(model, [IN(-1.0, 1.0) for _ in 1:2])
    x = [SV(vi) for vi in x_v]
    MOI.add_constraint(model, 1.0 * x[1] + 1.0 * x[2], LT(0.5))
    MOI.add_constraint(model, 1.0 * x[1] - 1.0 * x[2], LT(0.6))
    return model, x
end

function _is_equal(u::SAF, v::SAF)
    u_t = sort(
        u.terms,
        lt = (x, y) -> x.variable_index.value < y.variable_index.value,
    )
    v_t = sort(
        v.terms,
        lt = (x, y) -> x.variable_index.value < y.variable_index.value,
    )
    return (u_t == v_t) && (u.constant == v.constant)
end

@testset "naive_big_m_formulation!" begin
    for (activity_method, expected_big_ms) in [
        DisjunctiveConstraints.IntervalArithmetic() => (2.0, 2.0, -2.0, -2.0),
        DisjunctiveConstraints.LinearProgrammingRelaxation(Cbc.Optimizer) =>
            (0.5, 0.6, -2.0, -2.0),
    ]
        method = DisjunctiveConstraints.NaiveBigM(activity_method)
        model, x = _build_base_model()

        # Disjunction: (x_1 + x_2 <= 0) OR
        #              (x_1 - x_2 <= 0) OR
        #              (x_1 + 0.5x_2 >= 0.5 AND x_1 - 0.5x_2 >= 0.5)
        f_1 = 1.0 * x[1] + 1.0 * x[2]
        f_2 = 1.0 * x[1] - 1.0 * x[2]
        f_3 = 1.0 * x[1] + 0.5 * x[2]
        f_4 = 1.0 * x[1] - 0.5 * x[2]
        f = MOIU.vectorize([f_1, f_2, f_3, f_4])

        lbs = [
            -Inf -Inf -Inf
            -Inf -Inf -Inf
            -Inf -Inf 0.5
            -Inf -Inf 0.5
        ]
        ubs = [
            0.0 Inf Inf
            Inf 0.0 Inf
            Inf Inf Inf
            Inf Inf Inf
        ]
        s = DisjunctiveConstraints.DisjunctiveSet(lbs, ubs)
        disjunction = DisjunctiveConstraints.Disjunction(f, s)

        sum_ci, lt_cis, gt_cis =
            @inferred DisjunctiveConstraints.naive_big_m_formulation!(
                model,
                method,
                disjunction,
            )
        @test size(lt_cis) == (4, 3)
        @test size(gt_cis) == (4, 3)

        @test MOI.get(model, MOI.NumberOfVariables()) == 5
        for i in 1:2
            @test MOIU.get_bounds(model, Float64, VI(i)) == (-1.0, 1.0)
        end
        for i in 3:5
            @test MOI.is_valid(model, MOI.ConstraintIndex{SV,ZO}(i))
        end
        z = [SV(VI(i)) for i in 3:5]
        moi_et_cis = MOI.get(model, MOI.ListOfConstraintIndices{SAF,ET}())
        @test length(moi_et_cis) == 1
        @test sum_ci == moi_et_cis[1]
        let f = MOI.get(model, MOI.ConstraintFunction(), sum_ci),
            s = MOI.get(model, MOI.ConstraintSet(), sum_ci)

            @test _is_equal(f, 1.0 * z[1] + 1.0 * z[2] + 1.0 * z[3])
            @test s == ET(1.0)
        end
        moi_lt_cis = MOI.get(model, MOI.ListOfConstraintIndices{SAF,LT}())
        # two from base model, two from big-M formulation
        @test length(moi_lt_cis) == 2 + 2
        @test count(!isnothing, lt_cis) == 2
        let lt_ci = lt_cis[1, 1]
            @test !isnothing(lt_ci)
            @test lt_ci in moi_lt_cis
            f = MOI.get(model, MOI.ConstraintFunction(), lt_ci)
            s = MOI.get(model, MOI.ConstraintSet(), lt_ci)
            @test _is_equal(
                f,
                1.0 * x[1] + 1.0 * x[2] - expected_big_ms[1] * z[1],
            )
            @test s == MOI.LessThan(0.0)
        end

        let lt_ci = lt_cis[2, 2]
            @test !isnothing(lt_ci)
            @test lt_ci in moi_lt_cis
            f = MOI.get(model, MOI.ConstraintFunction(), lt_ci)
            s = MOI.get(model, MOI.ConstraintSet(), lt_ci)
            @test _is_equal(
                f,
                1.0 * x[1] - 1.0 * x[2] - expected_big_ms[2] * z[2],
            )
            @test s == MOI.LessThan(0.0)
        end

        moi_lt_cis = MOI.get(model, MOI.ListOfConstraintIndices{SAF,GT}())
        # None from base model, two from big-M formulation
        @test length(moi_lt_cis) == 2
        @test count(!isnothing, gt_cis) == 2
        let gt_ci = gt_cis[3, 3]
            @test !isnothing(gt_ci)
            @test gt_ci in moi_lt_cis
            f = MOI.get(model, MOI.ConstraintFunction(), gt_ci)
            s = MOI.get(model, MOI.ConstraintSet(), gt_ci)
            @test _is_equal(
                f,
                1.0 * x[1] + 0.5 * x[2] - expected_big_ms[3] * z[3],
            )
            @test s == MOI.GreaterThan(0.5)
        end

        let gt_ci = gt_cis[4, 3]
            @test !isnothing(gt_ci)
            @test gt_ci in moi_lt_cis
            f = MOI.get(model, MOI.ConstraintFunction(), gt_ci)
            s = MOI.get(model, MOI.ConstraintSet(), gt_ci)
            @test _is_equal(
                f,
                1.0 * x[1] - 0.5 * x[2] - expected_big_ms[3] * z[3],
            )
            @test s == MOI.GreaterThan(0.5)
        end
    end
end
