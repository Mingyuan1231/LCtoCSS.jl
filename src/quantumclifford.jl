"""
    CSSResult

Signed conversion returned by [`lc_to_css`](@ref).

- `checks`: pure-X rows followed by pure-Z rows, with physical signs.
- `hx`, `hz`: binary check matrices; `size(hx,1)` splits `checks` by type.
- `circuit`: QuantumClifford single-qubit operations in execution order.
- `transformed`: the original generators conjugated by the entire circuit.
- `generator_map`: each check is the signed product of the indicated rows of
  `transformed`. Thus exact syndromes transform as `s_css = K*s_input` over GF(2).
- `witness`: the binary local maps, input rank, factors, and anchor branches.
- `correction`: the final local Pauli correction, identity unless `positive=true`.

The circuit changes the local frame; the generator map changes the generating
set. Applying the circuit alone need not produce individually pure rows.
"""
struct CSSResult{S<:QC.Stabilizer,P<:QC.PauliOperator}
    checks::S
    hx::BitMatrix
    hz::BitMatrix
    circuit::Vector{QC.SingleQubitOperator}
    transformed::S
    generator_map::SparseMatrixCSC{Bool,Int}
    witness::CSSWitness
    correction::P
end

function _stabilizer(phases, x, z)
    if isempty(x)
        s = zero(QC.Stabilizer, size(x, 1), size(x, 2))
        QC.phases(s) .= phases
        return s
    end
    QC.Stabilizer(UInt8.(phases), x, z)
end

function _prepare(s::QC.Stabilizer)
    a = QC.stab_to_gf2(s)
    n = QC.nqubits(s)
    x, z = @view(a[:, 1:n]), @view(a[:, n+1:2n])
    g = _rref(_interleave(x, z))[1]
    _check_isotropic(g)
    all(p -> p == 0x00 || p == 0x02, QC.phases(s)) ||
        throw(ArgumentError("Stabilizer generators must have real phases (+1 or -1)"))
    # Own the storage even when the input is a view or uses another CPU layout.
    owned = _stabilizer(QC.phases(s), x, z)
    reduced, _, rank = QC.canonicalize!(copy(owned); ranks=true)
    any(!iszero, @view(QC.phases(reduced)[rank+1:end])) &&
        throw(ArgumentError("Inconsistent stabilizer generators: the group contains -I"))
    owned, g
end

# These words are in execution order. For example (:S,:H) implements U=H*S.
const _GATE_WORDS = Dict{LocalMap,Tuple{Vararg{Symbol}}}(
    (1, 0, 0, 1) => (),
    (0, 1, 1, 0) => (:H,),
    (1, 1, 0, 1) => (:S,),
    (1, 1, 1, 0) => (:S, :H),
    (0, 1, 1, 1) => (:H, :S),
    (1, 0, 1, 1) => (:H, :S, :H),
)

function _circuit(maps)
    circuit = QC.SingleQubitOperator[]
    for (i, q) in enumerate(maps), gate in _GATE_WORDS[q]
        push!(circuit, QC.SingleQubitOperator(gate == :H ? QC.sHadamard(i) : QC.sPhase(i)))
    end
    circuit
end

function _apply_circuit!(s, circuit)
    # QuantumClifford's packed kernels do not need to be called on empty tableaux.
    isempty(s) && return s
    for gate in circuit
        QC.apply!(s, gate)
    end
    s
end

function _projected_checks(t, independent)
    a = QC.stab_to_gf2(t)
    m, n = size(t)
    x, z = BitMatrix(a[:, 1:n]), BitMatrix(a[:, n+1:2n])
    ix = [i for i in 1:m if any(@view x[i, :])]
    iz = [i for i in 1:m if any(@view z[i, :])]
    if independent
        ix = ix[_rref(x[ix, :])[3]]
        iz = iz[_rref(z[iz, :])[3]]
    end
    hx, hz = x[ix, :], z[iz, :]
    checks = _stabilizer(zeros(UInt8, length(ix) + length(iz)),
                        vcat(hx, falses(length(iz), n)),
                        vcat(falses(length(ix), n), hz))
    hx, hz, checks
end

function _recover_products!(checks, transformed)
    a = QC.stab_to_gf2(transformed)
    selected = _rref(a)[3]
    g, pivots, _, coefficients = _rref(a[selected, :]; track=true)
    c = QC.stab_to_gf2(checks)
    row_indices, column_indices = Int[], Int[]
    for i in axes(c, 1)
        remaining = copy(c[i, :])
        combination = falses(length(selected))
        for (j, pivot) in enumerate(pivots)
            if remaining[pivot]
                remaining .⊻= @view g[j, :]
                combination .⊻= @view coefficients[j, :]
            end
        end
        any(remaining) && error("Internal error: a projected check is outside the transformed space")
        product = zero(QC.PauliOperator, QC.nqubits(transformed))
        for j in eachindex(selected)
            if combination[j]
                original = selected[j]
                push!(row_indices, i)
                push!(column_indices, original)
                QC.mul_left!(product, transformed[original])
            end
        end
        checks[i] = product
    end
    sparse(row_indices, column_indices, trues(length(row_indices)), size(c, 1), size(a, 1))
end

"""
    lc_to_css(s::QuantumClifford.Stabilizer; positive=false, independent=true)

Return a [`CSSResult`](@ref), or `nothing` when no local Clifford can make
the valid input CSS. Never mutate `s`. Invalid noncommuting, non-Hermitian,
or inconsistent signed generators raise `ArgumentError` (including a
redundant relation that generates `-I`).

By default return an independent subset of the nonzero X/Z projections of
the original transformed checks. This preserves the sparse-check weight
bound. With `independent=false`, retain all nonzero projections, including
duplicates: at most twice as many checks and at most twice the qubit incidence.

Signs are recovered by actual QuantumClifford Pauli multiplication. With
`positive=true`, append a product of single-qubit Paulis making every output
check positive. `result.circuit` includes that correction. For mixed or
destabilizer tableaux, pass `stabilizerview(s)` explicitly; logical operators
can then be conjugated by the returned circuit.

This constructs a CSS representative with deterministic tie-breaking, not
a unique invariant representative of an LC orbit.

# Example
```julia
using QuantumClifford, LCToCSS
s = S"XZ
      ZX"
result = lc_to_css(s; positive=true)
@assert result !== nothing
@assert size(result.hx, 1) + size(result.hz, 1) == 2
```
"""
function lc_to_css(s::QC.Stabilizer; positive::Bool=false, independent::Bool=true)
    transformed, g = _prepare(s)
    witness = _recognize(g)
    isnothing(witness) && return nothing
    circuit = _circuit(witness.local_maps)
    _apply_circuit!(transformed, circuit)
    hx, hz, checks = _projected_checks(transformed, independent)
    k = _recover_products!(checks, transformed)
    n = QC.nqubits(s)
    correction = zero(QC.PauliOperator, n)
    if positive
        c = QC.stab_to_gf2(checks)
        signs = BitVector(QC.phases(checks) .== 0x02)
        # C*Omega*t = signs in grouped (x | z) coordinates.
        bits = _solve(hcat(c[:, n+1:2n], c[:, 1:n]), signs)
        isnothing(bits) && error("Internal error: CSS phase correction is inconsistent")
        correction = QC.PauliOperator(bits[1:n], bits[n+1:2n])
        for i in 1:n
            x, z = bits[i], bits[n+i]
            !x && !z && continue
            gate = QC.SingleQubitOperator(x ? (z ? QC.sY(i) : QC.sX(i)) : QC.sZ(i))
            push!(circuit, gate)
            _apply_circuit!(checks, (gate,))
            _apply_circuit!(transformed, (gate,))
        end
    end
    CSSResult(checks, hx, hz, circuit, transformed, k, witness, correction)
end

function is_lc_css(s::QC.Stabilizer)
    _, g = _prepare(s)
    !isnothing(_recognize(g))
end
