module LCToCSS

import QuantumClifford as QC
using SparseArrays: SparseMatrixCSC, sparse

export CSSWitness, CSSResult, css_witness, is_lc_css, lc_to_css

include("gf2.jl")
include("recognition.jl")
include("quantumclifford.jl")

end
