function symplectic_int(a, b, n)
    isodd(sum(((a >> (2i)) & 1) * ((b >> (2i+1)) & 1) +
              ((b >> (2i)) & 1) * ((a >> (2i+1)) & 1) for i in 0:n-1; init=0))
end

function isotropic_spaces(n)
    spaces = Dict{Tuple{Vararg{Int}},Vector{Int}}((0,) => Int[])
    queue = Tuple{Vararg{Int}}[(0,)]
    cursor = 1
    while cursor <= length(queue)
        space = queue[cursor]
        cursor += 1
        generators = spaces[space]
        length(generators) == n && continue
        for v in 1:4^n-1
            v in space && continue
            any(g -> symplectic_int(g, v, n), generators) && continue
            expanded = Tuple(sort!(vcat(collect(space), [w ⊻ v for w in space])))
            haskey(spaces, expanded) && continue
            spaces[expanded] = vcat(generators, v)
            push!(queue, expanded)
        end
    end
    values(spaces)
end

# Literal search over all six invertible local maps. No projectors or anchors.
function literal_lc_oracle(x, z)
    n = size(x, 2)
    rank = reference_rank(hcat(x, z))
    maps = [(a, b, c, d) for a in (false, true), b in (false, true),
            c in (false, true), d in (false, true) if (a & d) ⊻ (b & c)]
    for choices in Iterators.product(ntuple(_ -> maps, n)...)
        tx, tz = binary_transform(x, z, choices)
        reference_rank(tx) + reference_rank(tz) == rank && return true
    end
    false
end

@testset "Literal LC search on every isotropic space through three qubits" begin
    total = 0
    for n in 0:3, rows in isotropic_spaces(n)
        x, z = unpack_rows(rows, n)
        witness = css_witness(x, z)
        @test !isnothing(witness) == literal_lc_oracle(x, z)
        isnothing(witness) || check_binary_witness(x, z, witness)
        total += 1
    end
    @test total == 550
end

function graph_matrix(mask, n)
    a = falses(n, n)
    k = 0
    for i in 1:n, j in i+1:n
        a[i, j] = a[j, i] = isodd(mask >> k)
        k += 1
    end
    a
end

function graph_mask(a)
    mask, k = 0, 0
    for i in axes(a, 1), j in i+1:size(a, 1)
        mask |= Int(a[i, j]) << k
        k += 1
    end
    mask
end

function bipartite(a)
    colors = fill(-1, size(a, 1))
    for root in eachindex(colors)
        colors[root] != -1 && continue
        colors[root] = 0
        queue = [root]
        while !isempty(queue)
            i = popfirst!(queue)
            for j in eachindex(colors)
                a[i, j] || continue
                colors[j] == colors[i] && return false
                if colors[j] == -1
                    colors[j] = 1 - colors[i]
                    push!(queue, j)
                end
            end
        end
    end
    true
end

@testset "Independent local-complementation orbits through five vertices" begin
    total = 0
    for n in 0:5
        count = 1 << (n*(n-1) ÷ 2)
        visited = falses(count)
        for start in 0:count-1
            visited[start+1] && continue
            orbit = Set([start])
            queue = [start]
            accepted = false
            while !isempty(queue)
                mask = popfirst!(queue)
                a = graph_matrix(mask, n)
                accepted |= bipartite(a)
                for v in 1:n
                    b = copy(a)
                    for i in 1:n, j in i+1:n
                        if a[v, i] && a[v, j]
                            b[i, j] = b[j, i] = !b[i, j]
                        end
                    end
                    next = graph_mask(b)
                    if !(next in orbit)
                        push!(orbit, next)
                        push!(queue, next)
                    end
                end
            end
            x = Matrix{Bool}(I, n, n)
            for mask in orbit
                visited[mask+1] = true
                z = graph_matrix(mask, n)
                w = css_witness(x, z)
                @test !isnothing(w) == accepted
                isnothing(w) || check_binary_witness(x, z, w)
                total += 1
            end
        end
        @test all(visited)
    end
    @test total == 1100
end

@testset "Second anchor branch and physical factors" begin
    # Truncated six-vertex path, with upper shear on vertex 3 moved to anchor 1.
    x = Matrix{Bool}(I, 6, 6)[1:5, :]
    z = [abs(i-j) == 1 for i in 1:5, j in 1:6]
    s = Stabilizer(x, z)
    apply!(s, sPhase(3))
    permutation = [3, 1, 2, 4, 5, 6]
    a = stab_to_gf2(s)
    s = Stabilizer(copy(QC.phases(s)), a[:, permutation], a[:, 6 .+ permutation])
    r = lc_to_css(s)
    @test r.witness.branches == [:c]
    verify_conversion(s, r)
    # Row mixing must not merge intrinsic factors.
    s = S"XXII\nZZII\nIIXX\nIIZZ"
    QC.mul_left!(s, 1, 3)
    QC.mul_left!(s, 2, 4)
    r = lc_to_css(s)
    @test r.witness.factors == [[1, 2], [3, 4]]
    verify_conversion(s, r)
end
