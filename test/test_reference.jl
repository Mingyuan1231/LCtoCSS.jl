@testset "Python reference decisions and witnesses" begin
    data = JSON3.read(read(joinpath(@__DIR__, "fixtures", "reference_cases.json"), String))
    @test length(data.cases) == 1136
    positives = 0
    for case in data.cases
        x, z = unpack_rows(case.rows, Int(case.n))
        w = css_witness(x, z)
        @test !isnothing(w) == case.equivalent
        if !isnothing(w)
            check_binary_witness(x, z, w)
            @test w.rank == case.rank
            positives += 1
        end
    end
    @test positives == 1034
end

@testset "Signed BB code integration: 72, 144, 288 qubits" begin
    data = JSON3.read(read(joinpath(@__DIR__, "fixtures", "bb_inputs.json"), String))
    for case in data.cases
        n = Int(case.n)
        x, z = unpack_rows(case.rows, n)
        s = Stabilizer(UInt8.(2 .* case.signs), x, z)
        before = copy(s)
        r = lc_to_css(s; positive=true, independent=false)
        verify_conversion(s, r; positive=true, independent=false)
        @test r.witness.rank == n - 12
        @test s == before
        @test maximum(sum(r.hx; dims=2); init=0) == 6
        @test maximum(sum(r.hz; dims=2); init=0) == 6
    end
end
