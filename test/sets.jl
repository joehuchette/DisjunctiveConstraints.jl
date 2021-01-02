@testset "DisjunctiveSet" begin
    for T in (Float64, Rational{Int})
        lbs = [T[1, 13//10], T[3//2, 4//3]]
        ubs = [T[5//2, 6], T[10, 11]]
        ds = DisjunctiveSet(lbs, ubs)
        @test ds.lbs == lbs
        @test ds.ubs == ubs
        @test MOI.dimension(ds) == 2

        bad_lbs = [T[1, 13//10, 2], T[3//2, 4//3]]
        @test_throws DimensionMismatch DisjunctiveSet(bad_lbs, ubs)
        bad_ubs = [T[5//2, 6], T[10, 11, 12]]
        @test_throws DimensionMismatch DisjunctiveSet(lbs, bad_ubs)
        @test_throws DimensionMismatch DisjunctiveSet(lbs, ubs[1:1])
    end
end

@testset "CombinatorialDisjunctiveSet" begin
    sparsity = [[1, 2], [2, 3], [3, 4]]
    cdc = CombinatorialDisjunctiveSet(sparsity)
    @test cdc.sparsity == sparsity
    @test MOI.dimension(cdc) == 4

    push!(sparsity, [1, 0])
    @test_throws ArgumentError CombinatorialDisjunctiveSet(sparsity)
end
