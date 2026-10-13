using QuantumClifford
using LCToCSS

# Generator changes alone cannot make XZ,ZX pure; one local H suffices.
input = S"XZ
          ZX"
result = lc_to_css(input; positive=true)
@assert result !== nothing
println("Input:")
println(input)
println("CSS checks:")
println(result.checks)
println("Local binary maps: ", result.witness.local_maps)
println("Generator map: ", result.generator_map)

# A rank-four code on five qubits with no CSS representative in its LC orbit.
five_qubit_code = S"XZZXI
                    IXZZX
                    XIXZZ
                    ZXIXZ"
@assert lc_to_css(five_qubit_code) === nothing
println("Five-qubit code rejected as expected.")

# Projection signs are determined by products, including dependent relations.
signed = S"YY
           ZZ
          -XX"
signed_result = lc_to_css(signed; independent=false)
@assert signed_result !== nothing
println("Signed CSS checks:")
println(signed_result.checks)
