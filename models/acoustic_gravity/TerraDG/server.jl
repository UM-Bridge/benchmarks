import Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.add(["UMBridge"])
using UMBridge
# TO EDIT: Path to TerraDG on your machine, script can be executed from anywhere
TerraDG_path="/home/areinarz/Desktop/TerraDG.jl"
Pkg.develop(path=TerraDG_path) # here an absolute path to the TerraDG installation, you can run the file from anywhere
using TerraDG

const CONFIG_FILE = TerraDG_path * "/src/earthquake.yaml"
const _base_config = TerraDG.Configuration(CONFIG_FILE)
const N_OUTPUT_STEPS = floor(Int, (_base_config.end_time - _base_config.pressure_sensors_start) / _base_config.pressure_sensors_step) + 1

function run_simulation(theta, sensor_locations)
    println("Julia threads available: ", Threads.nthreads())
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

    # Run simulation, with explicit unique paths
    t1 = time()
    output_prefix = joinpath(workdir, "output", "plot")
    overrides = Dict(
        "level_set_filename" => levelset_path,
        "output" => Dict(
            "pressure_sensors" => Dict(
                "locations" => sensor_locations
            )
        )
    )
    TerraDG.main(
        CONFIG_FILE;
        output_prefix=output_prefix,
        config_overrides=overrides,
    )
    print("Simulation done in ", time() - t1, ".\n Write output...")

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


function get_sensor_locations(config)
    raw_locations = config["captors"]
    return [(Float64(l[1]), Float64(l[2])) for l in raw_locations]
end

model = UMBridge.Model(
    name = "forward",
    inputSizes = [_base_config.grid_elements[1]],
    outputSizes = (config) -> [length(get_sensor_locations(config)) * N_OUTPUT_STEPS],
    evaluate = (input, config) -> [run_simulation(input[1], get_sensor_locations(config))]
)

UMBridge.serve_models([model], 4242)
