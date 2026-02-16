import Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.add(["UMBridge"])
using UMBridge
Pkg.develop(path="/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl") # here an absolute path to the TerraDG installation, you can run the file from anywhere
using TerraDG

function run_simulation(theta)
    # Write levelset input to file
    open("levelset.csv", write=true) do f
        write(f, "value\n") 
        for val in theta
            val = Int(val)
            write(f, "$(val)\n")
        end
    end

    # Run simulation
    TerraDG.main("/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl/src/earthquake.yaml")

    # Read pressure sensor output
    output = Float64[]
    open("/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl/output/plot_pressure_sensors.csv", read=true) do f
        readline(f) # skip header
        for line in eachline(f)
            parts = split(line, ",")
            push!(output, parse(Float64, parts[4]))
        end
    end
    print(output)
    return output
end

model = UMBridge.Model(
    name = "forward",
    inputSizes = [100],
    outputSizes = [101],
    evaluate = (input, config) -> [run_simulation(input[1])]
)

UMBridge.serve_models([model], 4242)
