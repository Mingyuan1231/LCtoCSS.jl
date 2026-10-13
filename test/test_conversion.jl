@testset "Boundary cases and invalid signed groups" begin
    for n in (0, 1, 5, 65), m in (0, 2), positive in (false, true)
        s = zero(Stabilizer, m, n)
        before = copy(s)
        r = lc_to_css(s; positive)
        verify_conversion(s, r; positive)
        @test s == before
        @test is_lc_css(s)
        @test all(==(LCToCSS.IDENTITY_MAP), r.witness.local_maps)
    end
    for s in (S"X\nZ", S"iX", S"-I", S"X\n-X", S"XX\nZZ\nYY")
        before = copy(s)
        @test_throws ArgumentError lc_to_css(s)
        @test_throws ArgumentError is_lc_css(s)
        @test s == before
    end
    @test_throws DimensionMismatch css_witness(falses(1, 2), falses(2, 1))
    @test_throws ArgumentError css_witness(Bool[1; 0;;], Bool[0; 1;;])
    bad = S"XZZXI\nIXZZX\nXIXZZ\nZXIXZ"
    @test !is_lc_css(bad)
    @test lc_to_css(bad) === nothing
    # An easy factor must not hide rejection in a second factor.
    product = S"XIIIII\nIXZZXI\nIIXZZX\nIXIXZZ\nIZXIXZ"
    @test lc_to_css(product) === nothing
    # A heavily redundant input should not allocate a dense m-by-m relation map.
    redundant = Stabilizer(trues(1000, 1), falses(1000, 1))
    result = lc_to_css(redundant; independent=false)
    @test nnz(result.generator_map) == 1000
    @test size(result.generator_map) == (1000, 1000)
end

@testset "Signed QuantumClifford conversion and sparse checks" begin
    rng = StableRNG(8713)
    families = [S"Y", S"YZ", S"XZ\nZX", S"YY\nZZ", S"XXX\nZZI\nIZZ", S"ZZI\nIZZ"]
    for base in families, shift in 0:5, positive in (false, true), independent in (false, true)
        s = copy(base)
        for i in eachindex(s)
            QC.phases(s)[i] = rand(rng, (0x00, 0x02))
        end
        maps = sort!(collect(keys(LCToCSS._GATE_WORDS)))
        for (i, q) in enumerate([maps[mod1(shift+j, 6)] for j in 1:nqubits(s)])
            for word in LCToCSS._GATE_WORDS[q]
                apply!(s, word == :H ? sHadamard(i) : sPhase(i))
            end
        end
        # Include a genuine signed redundant product and +I.
        if length(s) >= 2
            p = s[1] * s[2]
            s = vcat(s, Stabilizer([p, zero(PauliOperator, nqubits(s))]))
        end
        before = copy(s)
        result = lc_to_css(s; positive, independent)
        verify_conversion(s, result; positive, independent, rng)
        @test s == before
    end
    # A dependent sign relation is meaningful: YY*ZZ = -XX.
    s = S"YY\nZZ\n-XX"
    verify_conversion(s, lc_to_css(s))
    # stabilizerview is accepted, without changing logical/destabilizer storage.
    mixed = MixedDestabilizer(S"ZZI\nIZZ")
    before = copy(mixed)
    verify_conversion(stabilizerview(mixed), lc_to_css(stabilizerview(mixed)))
    @test mixed == before
end

@testset "All six lifts against physical 2x2 unitary matrices" begin
    id = Matrix{ComplexF64}(I, 2, 2)
    h = ComplexF64[1 1; 1 -1] / sqrt(2)
    phase = ComplexF64[1 0; 0 im]
    x, z = ComplexF64[0 1; 1 0], ComplexF64[1 0; 0 -1]
    paulis = Dict((false, false) => id, (true, false) => x,
                  (false, true) => z, (true, true) => im*x*z)
    for (q, words) in LCToCSS._GATE_WORDS
        u = copy(id)
        for word in words
            u = (word == :H ? h : phase) * u
        end
        for ((xb, zb), p) in paulis, sign in (0x00, 0x02)
            s = Stabilizer([sign], fill(xb, 1, 1), fill(zb, 1, 1))
            for gate in LCToCSS._circuit([q])
                apply!(s, gate)
            end
            actual = (QC.phases(s)[1] == 0 ? 1 : -1) * paulis[s[1, 1]]
            @test actual ≈ (sign == 0 ? 1 : -1) * u*p*u'
            @test s[1, 1] == ((xb & q[1]) ⊻ (zb & q[3]), (xb & q[2]) ⊻ (zb & q[4]))
        end
    end
end

@testset "Connected inputs across packed-word boundaries" begin
    for n in (31, 32, 33, 63, 64, 65, 128)
        x, z = falses(n, n), falses(n, n)
        x[1, :] .= true
        for i in 2:n
            z[i, i-1] = z[i, i] = true
        end
        s = Stabilizer(x, z)
        for i in 1:n
            apply!(s, sPhase(i))
            isodd(i) && apply!(s, sHadamard(i))
        end
        r = lc_to_css(s; positive=true)
        verify_conversion(s, r; positive=true)
        @test length(r.witness.factors) == 1
    end
end
