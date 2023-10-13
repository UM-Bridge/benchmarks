import umbridge
import json
import os
import csv

class AchlysModel(umbridge.Model):
    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [5] # Needs changing

    def get_output_sizes(self, config):
        return [500] # Needs changing

    def __call__(self, parameters, config):

        input = f"""
            {{
                &kt_grids_knobs
                    grid_option = 'single'
                /

                &kt_grids_single_parameters
                    aky = 0.7852773972394589
                    theta0 = 0.0
                /

                &theta_grid_parameters
                    ntheta = 128
                    nperiod = 9
                    akappa = 2.66
                    akappri = -0.25
                    rhoc = 0.67
                    tri = 0.3400000000000001
                    tripri = 0.25
                    shat = 0.87
                    qinp = 8.903370045724728
                    r_geo = 1.84
                    rmaj = 1.84
                    shift = -0.44
                    geotype = 0
                /

                &theta_grid_knobs
                    equilibrium_option = 'eik'
                /

                &theta_grid_eik_knobs
                    irho = 2
                    iflux = 0
                    s_hat_input = 3.1173933836235315
                    beta_prime_input = -0.64
                    local_eq = .true.
                    bishop = 4
                    equal_arc = .false.
                    writelots = .true.
                /

                &le_grids_knobs
                    ngauss = 8
                    negrid = 16
                    bouncefuzz = 1e-08
                /

                &dist_fn_knobs
                    adiabatic_option = 'iphi00=2'
                    opt_source = .true.
                    def_parity = .true.
                    even = .false.
                /

                &fields_knobs
                    field_option = 'implicit'
                    response_dir = 'response'
                /

                &knobs
                    fphi = 1.0
                    fapar = 1.0
                    fbpar = 0.0
                    delt = 0.02
                    nstep = 80000
                    wstar_units = .true.
                /

                &layouts_knobs
                    layout = 'xyles'
                /

                &collisions_knobs
                    collision_model = 'default'
                /

                &species_knobs
                    nspec = 3
                /

                &species_parameters_1
                    z = 1.0
                    mass = 1.0
                    dens = 0.5
                    temp = 1.0593220338983051
                    tprim = 0.0
                    fprim = 0.1450502618121151
                    uprim = 0.0
                    vnewk = 0.00039999999999999996
                    type = 'ion'
                    bess_fac = 1.0
                /

                &dist_fn_species_knobs_1
                    fexpr = 0.48
                    bakdif = 0.05
                /

                &species_parameters_2
                    z = 1.0
                    mass = 1.4975
                    dens = 0.5
                    temp = 1.0593220338983051
                    tprim = 0.0
                    fprim = 0.1450502618121151
                    uprim = 0.0
                    vnewk = 0.0003268711385781747
                    type = 'ion'
                    bess_fac = 1.0
                /

                &dist_fn_species_knobs_2
                    fexpr = 0.48
                    bakdif = 0.05
                /

                &species_parameters_3
                    z = -1
                    mass = 0.0002724
                    dens = 1.0
                    temp = 1.0
                    tprim = 4.98536245337014
                    fprim = 2.0329154559480602
                    uprim = 0.0
                    vnewk = 0.011195269013375201
                    type = 'electron'
                    bess_fac = 1.0
                /

                &dist_fn_species_knobs_3
                    fexpr = 0.48
                    bakdif = 0.05
                /

                &init_g_knobs
                    ginit_option = 'default_odd'
                    chop_side = .false.
                    phiinit = 1e-05
                    restart_dir = 'restart'
                /

                &gs2_diagnostics_knobs
                    write_ascii = .true.
                    write_omega = .true.
                    write_final_fields = .true.
                    write_fields = .true.
                    write_phi_over_time = .true.
                    write_bpar_over_time = .true.
                    write_apar_over_time = .true.
                    write_nl_flux_dist = .true.
                    write_fluxes = .true.
                    write_nl_flux = .true.
                    nwrite = 30
                    navg = 10
                    omegatol = 0.0001
                    omegatinst = 500.0
                    nsave = 5000
                    save_for_restart = .true.
                    write_final_epar = .true.
                /

                &parameters
                    beta = 0.29798123734456977
                    tite = 1.0
                    zeff = 1.0
                /

                &diagnostics_config
                    nwrite = 100000000
                /
            }}
        """

        # Write input to input file
        with open("/opt/fast.in", "w") as f:
            f.write(input)

        # add MOOSE Python libraries to PYTHONPATH, modify input files and run achlys
        os.system("sbatch some_script.sh") # submit to slurm or run directly?
        
        # Read results CSV file, write rows to output
        output = []
        with open("/opt/achlys/problems/thermal_desorption/ogorodnikova/tds_multiapp/desorp_multi_out.csv", "r") as f:
            reader = csv.reader(f)
            for row in reader:
                output.append(row)
                print(f"Read CSV Row: {row}")

        # Probably no need for interpolation
        # Interpolate data to a grid of size 500 before returning
        from scipy.interpolate import interp1d
        import pandas as pd
        from numpy import linspace, ndarray

        result = pd.read_csv('/opt/achlys/problems/thermal_desorption/ogorodnikova/tds_multiapp/desorp_multi_out.csv',usecols=['time','pfc_flux'])
        f = interp1d(result['time'],result['pfc_flux'])
        output = [ndarray.tolist(f(linspace(0,62.5,500)))]

        print(f"output: {output}")

        return output

    def supports_evaluate(self):
        return True

model = AchlysModel()

umbridge.serve_models([model], 4242)
