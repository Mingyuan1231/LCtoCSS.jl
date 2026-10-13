# QuantumClifford interface

## Conversion

`lc_to_css(stabilizer)` finds a local Clifford frame in which a valid stabilizer
group admits pure-X and pure-Z generators. It handles nonmaximal codes as well
as states. A result includes the physical circuit, signed CSS checks, and the
generator map needed for syndrome conversion.

This is distinct from [tableau canonicalization](https://qc.quantumsavory.org/stable/canonicalization/):
row operations expose structure in a fixed local frame. The conversion chooses
that frame by solving the factorwise affine preservation systems. The result
is a deterministic representative, with no assertion of uniqueness across an
LC orbit. `lc_to_css` therefore describes the guarantee more precisely than a
name implying a complete LC canonical form.

## Source boundaries

| File | Responsibility | QuantumClifford dependency |
| --- | --- | --- |
| `src/gf2.jl` | Packed ordinary binary elimination, affine solving, coordinates | None |
| `src/recognition.jl` | Physical factors, both affine branches, binary witness | None |
| `src/quantumclifford.jl` | Signed validation, Clifford lifts, sparse extraction, generator products and correction | Native tableaux, gates and Pauli multiplication |
| `src/LCToCSS.jl` | Package boundary and exports | Imports QuantumClifford |

The adapter uses `stab_to_gf2`, `canonicalize!`, `phases`, `Stabilizer`,
`PauliOperator`, `SingleQubitOperator`, `sHadamard`, `sPhase`, `sX`, `sY`, `sZ`,
`apply!`, and `mul_left!`. It owns its working tableau and does not depend on
QuantumClifford's internal bit-buffer layout. `canonicalize!` and `mul_left!`
handle genuine signed stabilizers only. Auxiliary affine systems use the
ordinary GF(2) solver.

## Return semantics

1. `nothing` is a mathematical rejection of a valid input; invalid stabilizers
   raise `ArgumentError`.
2. Circuit application and generator recombination are separate operations.
3. `positive=true` includes the Pauli correction in the returned circuit.
4. The generator map refers to the supplied input row order, including redundant
   rows, and acts on exact signed syndromes over GF(2).
5. No physical qubit permutation occurs. A view of the stabilizer portion can
   be supplied without modifying its surrounding mixed/destabilizer tableau.

The output row count can differ from the input. The nonmutating entry point
returns a separate stabilizer tableau and a circuit that callers can apply
to logical operators or other tableaux in the same frame.

LCToCSS defines its own functions and result types. It does not add methods
to QuantumClifford-owned functions. The binary helpers are private to the
LCToCSS module; the native tableau operations are confined to the adapter.

## Scope

The public binary interface takes two Boolean matrices. Packed words are an
internal representation. Sparse generators are preserved by selecting from
their projections; a sparse input storage type is not required. The package
executes on CPU and does not claim GPU support.

The general recognizer also accepts graph-state tableaux. It returns signed
CSS checks and a local Clifford circuit. Adjacency-specialized recognition,
explicit bipartite graph output, and serialized negative certificates are
outside this interface.
