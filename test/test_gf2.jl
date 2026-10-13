@testset "Packed GF(2) algebra and word boundaries" begin
    rng = StableRNG(304)
    for width in (0, 1, 2, 63, 64, 65, 127, 128, 129), rows in (0, 1, 8, 70)
        a = rand(rng, Bool, rows, width)
        g, pivots, selected, k = LCToCSS._rref(a; track=true)
        @test g == reference_rref(a)
        @test g == mod2(Int.(k) * Int.(a))
        @test length(pivots) == length(selected) == reference_rank(a[selected, :])
        @test reference_rref(a[selected, :]) == g
        @test LCToCSS._unpack(LCToCSS._pack(a), width) == a
        solution = rand(rng, Bool, width)
        rhs = mod2(Int.(a) * Int.(solution))
        answer = LCToCSS._solve(a, rhs)
        @test answer !== nothing
        @test mod2(Int.(a) * Int.(answer)) == rhs
        @test LCToCSS._solve(vcat(a, falses(1, width)), vcat(rhs, true)) === nothing
    end
    # These are ordinary linear equations, not an isotropic tableau.
    @test LCToCSS._solve(Bool[1 0; 0 1], Bool[1, 0]) == Bool[1, 0]
end
