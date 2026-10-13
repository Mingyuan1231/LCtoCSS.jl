# Testing

Run the test suite and example from the package root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. examples/basic.jl
```

The tests use local fixtures and seeded random inputs. They require neither
Python nor access to the companion repository.

## Coverage

| Test file | Coverage |
| --- | --- |
| `test_gf2.jl` | Packed GF(2) operations against scalar reference elimination |
| `test_conversion.jl` | Signed groups, native gates, generator maps, and input validation |
| `test_oracles.jl` | Independent exhaustive searches, anchor branches, and physical factors |
| `test_reference.jl` | Saved reference decisions and signed BB code inputs |

The suite checks:

- Packed elimination against separate scalar reference elimination, including
  tracked row relations, consistent and inconsistent systems, and machine-word
  boundaries.
- Literal six-choice LC search over all 550 isotropic spaces through three
  qubits, without using the affine criterion.
- Independent local-complementation orbit exploration over all 1,100 labeled
  graphs through five vertices, comparing orbit bipartiteness with recognition.
- All 1,136 reference cases, including independent verification of the 1,034
  reconstructed positive witnesses.
- Exact signed QuantumClifford group equality, generator reconstruction,
  syndrome transport, support bounds, zero spaces, invalid groups containing
  `-I`, both anchor branches, independent physical factors, and mixed-tableau
  stabilizer views.
- All six gate lifts against physical 2x2 complex unitary matrices.
- Connected instances across qubit word boundaries through 128 qubits, a
  1,000-row redundant input with a sparse generator map, and signed BB inputs
  on 72, 144 and 288 qubits.

The [fixture provenance](../test/fixtures/README.md) describes source hashes,
coordinate conventions, and the import command. Reference fixtures contain
input labels and expected decisions; the tests reconstruct and verify their
own witnesses.

## Dependency compatibility

To test the minimum supported QuantumClifford version in a local checkout:

```julia
using Pkg
Pkg.add(PackageSpec(name="QuantumClifford", version="0.11.4"))
Pkg.test()
```

The finite tests check implementation behavior and compatibility. The
polynomial-time correctness argument is the projector and factor theorem
in the accompanying paper.
