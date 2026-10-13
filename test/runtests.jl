using Test
using Random
using LinearAlgebra
import LinearAlgebra: I
using SparseArrays
using StableRNGs
using JSON3
using QuantumClifford
using LCToCSS

const QC = QuantumClifford
include("helpers.jl")

@testset "LCToCSS" begin
    include("test_gf2.jl")
    include("test_conversion.jl")
    include("test_oracles.jl")
    include("test_reference.jl")
end
