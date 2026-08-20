import Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.add(["UMBridge"])
using UMBridge
TerraDG_path="/home/nxdj93/TerraDG.jl"
Pkg.develop(path=TerraDG_path) # here an absolute path to the TerraDG installation, you can run the file from anywhere
using TerraDG

const CONFIG_FILE = TerraDG_path * "/src/earthquake.yaml"

function run_simulation(theta)
    print("Julia threads available: ", Threads.nthreads())
    workdir = mktempdir()

    t0 = time()
    print("Reading input...")
    levelset_path = joinpath(workdir, "levelset.csv")
    open(levelset_path, write=true) do f
        write(f, "value\n")
        for val in theta
            write(f, "$(Int(val))\n")
        end
    end
    elapsed = time() - t0
    print("Reading done in ", elapsed, ".\n Start simulation...")

    # Run simulation, with explicit unique paths - no reliance on cwd,
    # no shared files with any other concurrent call.
    t1 = time()
    output_prefix = joinpath(workdir, "output", "plot")
    TerraDG.main(
		         CONFIG_FILE;
			         output_prefix=output_prefix,
				         level_set_filename=levelset_path,
					     )
    elapsed = time() - t1
    print("Simulation done in ", elapsed, ".\n Write output...")

    # Read pressure sensor output
    t2 = time()
    output = Float64[]
    open(output_prefix * "_pressure_sensors.csv", read=true) do f
        readline(f) # skip header
        for line in eachline(f)
            parts = split(line, ",")
            push!(output, parse(Float64, parts[4]))
        end
    end
    elapsed = time() - t2
    print("Writing done in ", elapsed, ".\n")

    rm(workdir, recursive=true, force=true) # cleanup

    return output
end

model = UMBridge.Model(
    name = "forward",
    inputSizes = [100],
    outputSizes = [202],
    evaluate = (input, config) -> [run_simulation(input[1])]
)

UMBridge.serve_models([model], 4242)
