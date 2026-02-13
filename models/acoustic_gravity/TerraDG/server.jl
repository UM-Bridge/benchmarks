import Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.add(["UMBridge"])
using UMBridge

function run_simulation(theta)
    # Write levelset input to file
    open("levelset.csv", write=true) do f
        for val in theta
            write(f, "$(val)\n")
        end
    end

    # Run simulation
    include("TerraDG.jl")
    TerraDG.main("src/earthquake.yaml")

    # Read pressure sensor output
    output = Float64[]
    open("output/plot_pressure_sensors.csv", read=true) do f
        readline(f) # skip header
        for line in eachline(f)
            parts = split(line, ",")
            push!(output, parse(Float64, parts[4]))
        end
    end
    return output
end

model = UMBridge.Model(
    name = "forward",
    inputSizes = [100],
    outputSizes = [400],
    evaluate = (input, config) -> [run_simulation(input[1])]
)

UMBridge.serve_models([model], 4242)
