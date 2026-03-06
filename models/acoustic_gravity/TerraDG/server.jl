import Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.add(["UMBridge"])
using UMBridge
Pkg.develop(path="/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl") # here an absolute path to the TerraDG installation, you can run the file from anywhere
using TerraDG

function run_simulation(theta)
    t0 = time()
    print("Reading input...")
    # Write levelset input to file
    open("levelset.csv", write=true) do f
        write(f, "value\n") 
        for val in theta
            val = Int(val)
            write(f, "$(val)\n")
        end
    end
    elapsed = time()-t0
    print("Reading done in ", elapsed, ".\n Start simulation...")

    # Run simulation
    t1 = time()
    TerraDG.main("/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl/src/earthquake.yaml")

    elapsed = time()-t1
    print("Simulation done in ", elapsed,".\n Write output...")

    # Read pressure sensor output
    t2 = time()
    output = Float64[]
    open("/home/dubois/Dokumente/Code/Inversion/TerraDG/TerraDG.jl/output/plot_pressure_sensors.csv", read=true) do f
        readline(f) # skip header
        for line in eachline(f)
            parts = split(line, ",")
            push!(output, parse(Float64, parts[4]))
        end
    end

    elapsed = time()-t2
    print("Writing done in ", elapsed,".\n")
    return output
end

model = UMBridge.Model(
    name = "forward",
    inputSizes = [100],
    outputSizes = [202],
    evaluate = (input, config) -> [run_simulation(input[1])]
)

UMBridge.serve_models([model], 4242)
