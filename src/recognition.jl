const LocalMap = NTuple{4,Bool}
const IDENTITY_MAP = (true, false, false, true)

"""
    CSSWitness

A binary LC-to-CSS witness. `local_maps[i] = (a,b,c,d)` represents the
right-acting matrix `[a b; c d]`: `(x,z) -> (x*a ⊻ z*c, x*b ⊻ z*d)`.
`rank` is the input stabilizer rank. `factors` lists the physical qubits in
each indecomposable factor; `branches[i]` is `:b` or `:c`, the successful
anchor branch on its first qubit. Indices are one-based.
"""
struct CSSWitness
    local_maps::Vector{LocalMap}
    rank::Int
    factors::Vector{Vector{Int}}
    branches::Vector{Symbol}
end

function _factors(g)
    n = size(g, 2) ÷ 2
    parent = collect(1:n)
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        i
    end
    for row in eachrow(g)
        first = findfirst(row)
        isnothing(first) && continue
        anchor = (first + 1) ÷ 2
        for j in eachindex(row)
            row[j] && (parent[root((j + 1) ÷ 2)] = root(anchor))
        end
    end
    components = Vector{Int}[]
    index = Dict{Int,Int}()
    for i in 1:n
        k = get!(index, root(i)) do
            push!(components, Int[])
            length(components)
        end
        push!(components[k], i)
    end
    components
end

# Nonpivot equations of GP = (GP)[:,pivots]G, with an extra anchor row last.
# Each GP coordinate contains at most three terms, so construction is O(N^3).
function _preserving_equations(g, pivots)
    r, width = size(g)
    n = width ÷ 2
    nonpivots = setdiff(1:width, pivots)
    equations = zeros(UInt64, cld(3n + 1, 64), r * length(nonpivots) + 1)
    function add_coordinate!(eq, row, column)
        q = (column + 1) ÷ 2
        x, z = g[row, 2q-1], g[row, 2q]
        if isodd(column)
            x && _flip!(eq, 3q - 2)
            z && _flip!(eq, 3q)
        else
            x && _flip!(eq, 3q - 1)
            z && _flip!(eq, 3q - 2)
            z && _flip!(eq, 3n + 1)
        end
    end
    k = 0
    for row in 1:r, column in nonpivots
        k += 1
        eq = @view equations[:, k]
        add_coordinate!(eq, row, column)
        for (beta, pivot) in enumerate(pivots)
            g[beta, column] && add_coordinate!(eq, row, pivot)
        end
    end
    equations
end

function _local_map(a, b, c)
    b && c && error("Internal error: the affine solution is not a projector")
    vectors = ((true, false), (false, true), (true, true))
    action(v) = ((a & v[1]) ⊻ (b & v[2]), (c & v[1]) ⊻ (!a & v[2]))
    u = only(v for v in vectors if action(v) == v)
    w = only(v for v in vectors if action(v) == (false, false))
    (u[1], w[1], u[2], w[2])
end

function _recognize(g)
    n = size(g, 2) ÷ 2
    maps = fill(IDENTITY_MAP, n)
    components = _factors(g)
    branches = Symbol[]
    for component in components
        columns = [j for i in component for j in (2i-1, 2i)]
        localg, pivots = _rref(g[:, columns])
        if isempty(pivots)
            push!(branches, :b)
            continue
        end
        equations = _preserving_equations(localg, pivots)
        anchor = @view equations[:, end]
        _flip!(anchor, 2) # b_1 = 0
        solution = _solve(equations, 3length(component))
        branch = :b
        if isnothing(solution)
            _flip!(anchor, 2)
            _flip!(anchor, 3) # c_1 = 0
            solution = _solve(equations, 3length(component))
            branch = :c
        end
        isnothing(solution) && return nothing
        push!(branches, branch)
        for (j, i) in enumerate(component)
            maps[i] = _local_map(solution[3j-2], solution[3j-1], solution[3j])
        end
    end
    CSSWitness(maps, size(g, 1), components, branches)
end

"""
    css_witness(x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool})

Decide whether the binary stabilizer row space `(x | z)` is locally
Clifford-equivalent to CSS. Return a [`CSSWitness`](@ref), or `nothing` if
no local Clifford exists. Noncommuting inputs raise `ArgumentError`.
Signs are not part of this binary interface. Inputs are never modified.

The deterministic factorwise algorithm solves at most two affine systems
per physical factor. It supports redundant/zero rows, ranks below the qubit
count, unused qubits, and zero qubits. The scalar-operation bound is
O(m*n^2 + n^4); equations and elimination rows are packed in UInt64 words.
The chosen representative need not be shared by LC-equivalent inputs.
"""
function css_witness(x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool})
    g = _rref(_interleave(x, z))[1]
    _check_isotropic(g)
    _recognize(g)
end

"""
    is_lc_css(x, z)
    is_lc_css(stabilizer)

Return whether a CSS representative exists under single-qubit Cliffords.
Invalid inputs raise an exception; `false` means a valid input was rejected.
"""
is_lc_css(x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool}) = !isnothing(css_witness(x, z))
