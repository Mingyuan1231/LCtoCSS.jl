# Each logical row is one column of machine words, so row XORs are contiguous.
# This is ordinary GF(2) algebra: equations need not be commuting Pauli rows.
@inline _bit(a, j) = !iszero(a[(j - 1) ÷ 64 + 1] & (UInt64(1) << ((j - 1) % 64)))
@inline function _flip!(a, j)
    a[(j - 1) ÷ 64 + 1] ⊻= UInt64(1) << ((j - 1) % 64)
    a
end

function _xor!(a, b)
    @inbounds @simd for k in eachindex(a, b)
        a[k] ⊻= b[k]
    end
    a
end

function _pack(a::AbstractMatrix{Bool})
    rows, width = size(a)
    p = zeros(UInt64, cld(width, 64), rows)
    for i in 1:rows, j in 1:width
        a[i, j] && _flip!(@view(p[:, i]), j)
    end
    p
end

function _unpack(p, width)
    a = falses(size(p, 2), width)
    for i in axes(a, 1), j in axes(a, 2)
        a[i, j] = _bit(@view(p[:, i]), j)
    end
    a
end

function _swap!(p, i, j)
    @inbounds for k in axes(p, 1)
        p[k, i], p[k, j] = p[k, j], p[k, i]
    end
end

# selected contains independent original row indices, even though p is reduced.
# Tracking is only used on an already independent basis, avoiding an m-by-m map.
function _rref(a::AbstractMatrix{Bool}; track::Bool=false)
    m, width = size(a)
    p = _pack(a)
    coefficients = zeros(UInt64, track ? cld(m, 64) : 0, m)
    if track
        for i in 1:m
            _flip!(@view(coefficients[:, i]), i)
        end
    end
    original = collect(1:m)
    pivots = Int[]
    for j in 1:width
        k = length(pivots) + 1
        found = findfirst(i -> _bit(@view(p[:, i]), j), k:m)
        isnothing(found) && continue
        i = k + found - 1
        _swap!(p, k, i)
        _swap!(coefficients, k, i)
        original[k], original[i] = original[i], original[k]
        for row in 1:m
            if row != k && _bit(@view(p[:, row]), j)
                _xor!(@view(p[:, row]), @view(p[:, k]))
                track && _xor!(@view(coefficients[:, row]), @view(coefficients[:, k]))
            end
        end
        push!(pivots, j)
    end
    r = length(pivots)
    (_unpack(@view(p[:, 1:r]), width), pivots, original[1:r],
     track ? _unpack(@view(coefficients[:, 1:r]), m) : falses(0, 0))
end

_rank(a::AbstractMatrix{Bool}) = length(_rref(a)[2])

# Solve a packed augmented system; the RHS is bit width+1. Free variables are 0.
function _solve(equations::Matrix{UInt64}, width::Int)
    pivots = zeros(Int, width)
    echelon = zeros(UInt64, size(equations, 1), width)
    rank = 0
    value = zeros(UInt64, size(equations, 1))
    for row in axes(equations, 2)
        copyto!(value, @view(equations[:, row]))
        inserted = false
        for j in 1:width
            _bit(value, j) || continue
            if pivots[j] == 0
                rank += 1
                pivots[j] = rank
                echelon[:, rank] = value
                inserted = true
                break
            end
            _xor!(value, @view(echelon[:, pivots[j]]))
        end
        !inserted && _bit(value, width + 1) && return nothing
    end
    solution = falses(width)
    for j in width:-1:1
        pivots[j] == 0 && continue
        row = @view echelon[:, pivots[j]]
        value = _bit(row, width + 1)
        for k in j+1:width
            value ⊻= _bit(row, k) & solution[k]
        end
        solution[j] = value
    end
    solution
end

function _solve(a::AbstractMatrix{Bool}, b::AbstractVector{Bool})
    size(a, 1) == length(b) || throw(DimensionMismatch("One RHS bit is required per equation"))
    _solve(_pack(hcat(a, b)), size(a, 2))
end

function _interleave(x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool})
    size(x) == size(z) || throw(DimensionMismatch("X and Z matrices must have the same size"))
    Base.require_one_based_indexing(x, z)
    a = falses(size(x, 1), 2size(x, 2))
    a[:, 1:2:end] = x
    a[:, 2:2:end] = z
    a
end

function _check_isotropic(g)
    x, z = _pack(@view(g[:, 1:2:end])), _pack(@view(g[:, 2:2:end]))
    for i in axes(g, 1), j in 1:i-1
        parity = false
        @inbounds for k in axes(x, 1)
            parity ⊻= isodd(count_ones((x[k, i] & z[k, j]) ⊻ (z[k, i] & x[k, j])))
        end
        parity && throw(ArgumentError("Stabilizer generators must commute"))
    end
    nothing
end
