# Run with --project=test and the path to a companion LCtoCSS checkout.
using JSON3
using SHA

function main(args)
    length(args) == 1 || throw(ArgumentError("Usage: julia --project=test tools/import_reference_fixtures.jl /path/to/LCtoCSS"))
    source = abspath(only(args))
    destination = normpath(joinpath(@__DIR__, "..", "test", "fixtures"))
    mkpath(destination)
    file = joinpath(source, "validation", "cases.json")
    data = JSON3.read(read(file, String))
    cases = [(name=c.name, n=c.n, rank=c.rank, rows=["0x" * string(r; base=16) for r in c.input_rows],
              equivalent=c.oracle.equivalent) for c in data.cases]
    @assert all(c.projector.equivalent == c.oracle.equivalent for c in data.cases)
    open(joinpath(destination, "reference_cases.json"), "w") do io
        JSON3.write(io, (source="validation/cases.json", sha256=bytes2hex(sha256(read(file))), cases=cases))
        println(io)
    end
    bb = []
    for n in (72, 144, 288)
        relative = "application_demo/results/bb_$(n)_instance.json"
        file = joinpath(source, split(relative, '/')...)
        data = JSON3.read(read(file, String))
        push!(bb, (source=relative, sha256=bytes2hex(sha256(read(file))), n=n,
                   rows=data.input_rows_hex, signs=data.input_signs))
    end
    open(joinpath(destination, "bb_inputs.json"), "w") do io
        JSON3.write(io, (cases=bb,))
        println(io)
    end
    println("Imported $(length(cases)) binary fixtures and $(length(bb)) signed BB inputs.")
end

main(ARGS)
