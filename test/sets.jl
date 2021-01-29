@testset "DisjunctiveSet" begin
    for T in (Float64, Rational{Int})
        lbs = T[1 13//10 4; 3//2 4//3 7]
        ubs = T[5//2 6 2; 10 11 12]
        ds = @inferred DisjunctiveSet(lbs, ubs)
        @test ds.lbs == lbs
        @test ds.ubs == ubs
        @test @inferred MOI.dimension(ds) == 2
        @test @inferred DisjunctiveConstraints.num_alternatives(ds) == 3

        bad_lbs = T[1 13//10; 3//2 4//3]
        @test_throws DimensionMismatch DisjunctiveSet(bad_lbs, ubs)
    end
end

@testset "CombinatorialDisjunctiveSet" begin
    sparsity = [[1, 2], [2, 3], [3, 4]]
    cdc = @inferred CombinatorialDisjunctiveSet(sparsity)
    @test cdc.sparsity == sparsity
    @test @inferred MOI.dimension(cdc) == 4

    push!(sparsity, [1, 0])
    @test_throws ArgumentError CombinatorialDisjunctiveSet(sparsity)
end

@testset "Disjunction" begin
    f = MOI.VectorOfVariables([MOI.VariableIndex(1), MOI.VariableIndex(3)])
    lbs = [-1 -1 -1; -2 -2 -2]
    ubs = [1 1 1; 2 2 2]
    s = DisjunctiveConstraints.DisjunctiveSet(lbs, ubs)
    disjunction = @inferred DisjunctiveConstraints.Disjunction(f, s)
    @test @inferred DisjunctiveConstraints.num_alternatives(disjunction) == 3

    lbs = vcat(lbs, [-3 -3 -3])
    ubs = vcat(ubs, [3 3 3])
    s_bad = DisjunctiveConstraints.DisjunctiveSet(lbs, ubs)
    @test_throws DimensionMismatch DisjunctiveConstraints.Disjunction(f, s_bad)
end
