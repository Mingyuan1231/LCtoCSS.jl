# Deliberately scalar, independent reference elimination for test assertions.
function reference_rref(a)
    b = Matrix{Bool}(a)
    r = 0
    for j in axes(b, 2)
        p = findfirst(i -> b[i, j], r+1:size(b, 1))
        isnothing(p) && continue
        p += r
        r += 1
        b[r, :], b[p, :] = copy(b[p, :]), copy(b[r, :])
        for i in axes(b, 1)
            if i != r && b[i, j]
                for k in axes(b, 2)
                    b[i, k] ⊻= b[r, k]
                end
            end
        end
    end
    b[1:r, :]
end

reference_rank(a) = size(reference_rref(a), 1)
mod2(a) = isodd.(a)

function unpack_rows(rows, n)
    x, z = falses(length(rows), n), falses(length(rows), n)
    for (i, encoded) in enumerate(rows)
        row = encoded isa AbstractString ? parse(BigInt, encoded) : BigInt(encoded)
        for j in 1:n
            x[i, j] = isodd(row >> (2j - 2))
            z[i, j] = isodd(row >> (2j - 1))
        end
    end
    x, z
end

function binary_transform(x, z, maps)
    tx, tz = copy(x), copy(z)
    for (j, (a, b, c, d)) in enumerate(maps), i in axes(x, 1)
        tx[i, j] = (x[i, j] & a) ⊻ (z[i, j] & c)
        tz[i, j] = (x[i, j] & b) ⊻ (z[i, j] & d)
    end
    tx, tz
end

function check_binary_witness(x, z, w)
    @test all(q -> (q[1] & q[4]) ⊻ (q[2] & q[3]), w.local_maps)
    tx, tz = binary_transform(x, z, w.local_maps)
    @test reference_rank(tx) + reference_rank(tz) == w.rank == reference_rank(hcat(x, z))
    @test all(iszero, mod2(Int.(tx) * Int.(tz)'))
end

function group_basis(s)
    a, _, r = canonicalize!(copy(s); ranks=true)
    a[1:r]
end

function verify_conversion(s, result; positive=false, independent=true, rng=StableRNG(1))
    @test result !== nothing
    isnothing(result) && return
    r = result
    n, m = nqubits(s), length(s)
    transformed = copy(s)
    for gate in r.circuit
        isempty(transformed) || apply!(transformed, gate)
    end
    @test transformed == r.transformed
    @test group_basis(transformed) == group_basis(r.checks)
    @test size(r.hx, 2) == size(r.hz, 2) == n
    @test reference_rank(r.hx) + reference_rank(r.hz) == r.witness.rank
    @test all(iszero, mod2(Int.(r.hx) * Int.(r.hz)'))
    @test size(r.generator_map) == (length(r.checks), m)
    @test r.generator_map isa SparseMatrixCSC{Bool,Int}
    c = stab_to_gf2(r.checks)
    @test c == vcat(hcat(r.hx, falses(size(r.hx))), hcat(falses(size(r.hz)), r.hz))
    @test c == mod2(Int.(r.generator_map) * Int.(stab_to_gf2(transformed)))
    if positive
        @test all(iszero, QC.phases(r.checks))
    end
    if independent
        @test length(r.checks) == r.witness.rank
    end
    for i in eachindex(r.checks)
        product = zero(PauliOperator, n)
        for j in 1:m
            r.generator_map[i, j] && QC.mul_left!(product, transformed[j])
        end
        @test product == r.checks[i]
    end
    input_support = [count(j -> s[i, j] != (false, false), 1:n) for i in 1:m]
    output_support = [count(j -> r.checks[i, j] != (false, false), 1:n) for i in eachindex(r.checks)]
    @test maximum(output_support; init=0) <= maximum(input_support; init=0)
    @test length(r.checks) <= 2m
    for j in 1:n
        @test count(i -> r.checks[i, j] != (false, false), eachindex(r.checks)) <=
              2count(i -> s[i, j] != (false, false), 1:m)
    end
    # Exact syndrome transport for errors in the original coordinate frame.
    b = stab_to_gf2(s)
    for _ in 1:3
        ex, ez = rand(rng, Bool, 1, n), rand(rng, Bool, 1, n)
        tx, tz = binary_transform(ex, ez, r.witness.local_maps)
        input_syndrome = mod2(Int.(b[:, 1:n]) * vec(Int.(ez)) + Int.(b[:, n+1:2n]) * vec(Int.(ex)))
        css_syndrome = mod2(Int.(c[:, 1:n]) * vec(Int.(tz)) + Int.(c[:, n+1:2n]) * vec(Int.(tx)))
        @test css_syndrome == mod2(Int.(r.generator_map) * Int.(input_syndrome))
    end
end
