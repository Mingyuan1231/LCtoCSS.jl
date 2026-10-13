# Fixture provenance

These compact fixtures were extracted from the companion Python
[LCtoCSS repository](https://github.com/Mingyuan1231/LCtoCSS).
Each source's SHA-256 digest is embedded in the JSON output.

- `reference_cases.json`: all 1,136 inputs from `validation/cases.json`, with
  decisions from its independent exhaustive local-Clifford oracle. Stored
  projector decisions agree with the oracle. Only input labels, dimensions,
  ranks, names, and expected decisions are copied; the Julia implementation
  constructs its own witnesses.
- `bb_inputs.json`: supplied labels and signs from the recorded 72-, 144-,
  and 288-qubit BB experiments in `application_demo/results/`. Generation and
  recovery gates are omitted so the Julia recognizer sees only the input.

Labels are hexadecimal integers with interleaved x,z bits, least-significant
qubit first. Sign bits 0 and 1 mean +1 and -1. The Julia importer converts the
record format explicitly; QuantumClifford itself uses grouped x,z columns.

Regenerate from the companion repository with Julia, from this package root:

```sh
julia --project=test tools/import_reference_fixtures.jl /path/to/LCtoCSS
```

The source checkout path is required. The importer uses the test environment's
JSON3 and SHA packages. Ordinary package use and tests
read only the fixtures in this directory. The fixtures derive from original
MIT-licensed LCtoCSS validation data; no third-party catalog is redistributed.
