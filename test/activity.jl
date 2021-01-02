const FACTORY = () -> begin
    model = Cbc.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    return model
end
const OPTIMIZER = FACTORY()
const MODEL = MOI.Utilities.CachingOptimizer(
    MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    OPTIMIZER,
)

# Feasible region: {x in R^3 :
#     x[1] + x[2] + x[3] <= 1
#     2x[1] - x[3] >= 0
#     0 <= x <= 1
#     x[3] in {0,1}
function _repopulate_model!(model)
    MOI.empty!(model)
    _x = [
        MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
        for i in 1:3
    ]
    x = [MOI.SingleVariable(_x[i][1]) for i in 1:3]
    MOI.add_constraint(
        model,
        1.0 * x[1] + 1.0 * x[2] + 1.0 * x[3],
        MOI.LessThan(1.0),
    )
    MOI.add_constraint(model, 2.0 * x[1] - 1.0 * x[3], MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x[3], MOI.ZeroOne())
    return x
end

@testset "maximum_activity" begin
    @testset "IntervalArithmetic" begin
        x = _repopulate_model!(MODEL)

        @test maximum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.IntervalArithmetic(),
        ) ≈ 4.0
        @test minimum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.IntervalArithmetic(),
        ) ≈ 0.0
    end

    @testset "LinearProgrammingRelaxation" begin
        x = _repopulate_model!(MODEL)

        @test maximum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.LinearProgrammingRelaxation(FACTORY),
        ) ≈ 2.0
        @test minimum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.LinearProgrammingRelaxation(FACTORY),
        ) ≈ 0.0
    end

    @testset "FullFormulation" begin
        x = _repopulate_model!(MODEL)

        @test maximum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.FullFormulation(FACTORY),
        ) ≈ 1.0
        @test minimum_activity(
            MODEL,
            1.0 * x[2] + 3.0 * x[3],
            DisjunctiveConstraints.FullFormulation(FACTORY),
        ) ≈ 0.0
    end
end
