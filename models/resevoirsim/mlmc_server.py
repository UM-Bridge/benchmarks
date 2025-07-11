import umbridge
import subprocess
import os
import pickle
import logging
import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO)

class PFLOTRANModel(umbridge.Model):
    def __init__(self):
        super().__init__("pflotran_simulation")
        self.simulation_script = "~/pflotran_ogs_1.8/src/pflotran/pflotran"
        self.cnt = 0 # Count simulations

    def get_input_sizes(self, config):
        return [15]  # 12 parameters for each layer + perm + 2 faults

    def get_output_sizes(self, config):
        return [1]

    # Catch output error if nan
    # check how many errors
    def __call__(self, parameters, config):
        try:
            # Prepare the input file for PFLOTRAN based on parameters
            iteration = str(config.get("iteration"))
            self.output_folder = "~/output/"

            self.level = int(config.get("l"))

            print(np.array(parameters))

            """
            Introduced UM_JOB_DIR env variables to redirect stdout and stderr
            to respective iteration folder.
            Problem: Will have race condition if multiple server on on node
            """
            
            # GP or pflotran
            if self.level == 0:
                gp_folder = self.output_folder + "gp/" + iteration + os.sep
                os.environ["UM_JOB_DIR"] = gp_folder # Problematic if multiple server in one node
                os.system(f"mkdir -p {gp_folder}")
                result = self.gp(np.array(parameters).reshape(1, -1)).tolist()
                logging.info("Evaluated GP surrogate")

                gp_params = parameters + result

                with open(f"{gp_folder}GP.pkl", "wb") as h:
                    pickle.dump(gp_params, h)

            elif self.level == 1:  
                # Run the PFLOTRAN simulation
                self.pflotran_folder = self.output_folder + "pflotran/coarse/" + iteration + os.sep
                self.input_folder = "./coarse_model/"
                self.pflotran_input_path = self.pflotran_folder + "Coarse_CCS_3DF.in"
                self.infile = self.pflotran_folder + "Coarse_CCS_3DF-mas.dat"
                self.outfile = self.pflotran_folder + "Coarse_CCS_3DF.csv"#+"_"+str(self.cnt)+".csv"
                self.external_file = "ccs_3df_geom.grdecl"

                os.system(f"mkdir -p {self.pflotran_folder}")
                os.system(f"cp -Tr {self.input_folder} {self.pflotran_folder}")
                self.prepare_input_file(parameters[0])
                logging.info("Input file prepared successfully.")
                self.run_simulation()
                logging.info("Simulation ran successfully.")

                # Process the output file to extract results
                result = self.process_output_file()
                logging.info("Output processed successfully.")

            elif self.level == 2:  
                # Run the PFLOTRAN simulation
                self.pflotran_folder = self.output_folder + "pflotran/fine/" + iteration + os.sep
                self.input_folder = "./fine_model/"
                self.pflotran_input_path = self.pflotran_folder + "FINE_CCS_3DF.in"
                self.infile = self.pflotran_folder + "FINE_CCS_3DF-mas.dat"
                self.outfile = self.pflotran_folder + "FINE_CCS_3DF.csv"#+"_"+str(self.cnt)+".csv"
                self.external_file = "ccs_3df_geom_refined.grdecl"
                
                os.system(f"mkdir -p {self.pflotran_folder}")
                os.system(f"cp -Tr {self.input_folder} {self.pflotran_folder}")
                self.prepare_input_file(parameters[0])
                logging.info("Input file prepared successfully.")
                self.run_simulation()
                logging.info("Simulation ran successfully.")

                # Process the output file to extract results
                result = self.process_output_file()
                logging.info("Output processed successfully.")
            
            else:
                raise Exception("Invalid level")

            print(result)
        except Exception as e:
            logging.error(f"Error during simulation: {e}")
            raise
        
        else:
            return result

    def gp(self, X):
        with open("./opengosim_gp.pkl", "rb") as h:
            GP = pickle.load(h)
        with open("./opengosim_xscaler.pkl", "rb") as h:
            xscaler = pickle.load(h)
        with open("./opengosim_yscaler.pkl", "rb") as h:
            yscaler = pickle.load(h)
        X = xscaler.transform(X)
        output = yscaler.inverse_transform(GP.predict(X).reshape(-1, 1)) / 1e8
        return output
        
    def process_output_file(self):
        # Converts a -mas.dat into a .csv file
        replacements = {' "':",",'"':'','  ':',',' -':',-'}

        with open(self.infile) as infile, open(self.outfile, 'w') as outfile:
            for line in infile:
                for src, target in replacements.items():
                    line = line.replace(src, target)
                outfile.write(line[1:])

        # Extract feature from .csv file
        try:
            df=pd.read_csv(self.outfile)
        except Exception as e:
            print(f"Error reading output csv: {e}")
        else:
            result = []
            result.append(float(df["Field fgit [m^3]"].iloc[-1]) / 1e8) # Total gas injected. Scaled by 1e8 for numerical precision

            # We don't care below for now
            """
            result.append(float(df["BPressure [bar] 30_30_1"].iloc[-1])) # Pressure should never exceed 200
            result.append(float(df["BPressure [bar] 7_7_1"].iloc[-1]))
            result.append(float(df["BPressure [bar] 27_11_1"].iloc[-1]))
            result.append(float(df["BSgas 7_1_1"].iloc[-1])) # Gas in NOGO area, should remain zero
            result.append(float(df["BSgas 7_2_1"].iloc[-1]))
            result.append(float(df["BSgas 7_3_1"].iloc[-1]))
            result.append(float(df["BSgas 7_4_1"].iloc[-1]))
            result.append(float(df["BSgas 7_5_1"].iloc[-1]))
            result.append(float(df["BSgas 7_6_1"].iloc[-1]))
            result.append(float(df["BSgas 7_7_1"].iloc[-1]))
            result.append(float(df["BSgas 1_7_1"].iloc[-1]))
            result.append(float(df["BSgas 2_7_1"].iloc[-1]))
            result.append(float(df["BSgas 3_7_1"].iloc[-1]))
            result.append(float(df["BSgas 4_7_1"].iloc[-1]))
            result.append(float(df["BSgas 5_7_1"].iloc[-1]))
            result.append(float(df["BSgas 6_7_1"].iloc[-1]))
            result.append(float(df["BSgas 7_7_1"].iloc[-1]))
            """
        return [result]


    def run_simulation(self):
        self.cnt = self.cnt+1
        try:
            os.system("mpirun -n 6 " + self.simulation_script + " -pflotranin "+  self.pflotran_input_path)

        except Exception as e:
            logging.error(f"Error running simulation: {e}")
            raise

    def prepare_input_file(self, parameters):
        try:
            porosities = parameters[:12]
            permz = parameters[12]
            faulta = parameters[13]
            faultb = parameters[14]
            l = self.level  # self.level is 0 indexed
            odd_sequence = lambda x: (2 ** (l - 1)) * (i - 1) + 1
            even_sequence = lambda x: (2 ** (l - 1)) * i
            input_data = f"""DIMENS
{30*l} {30*l} {12*l} /

external_file {self.external_file}

-- Control zone
EQUALS
FIPNOGO 5 1 {7*l} 1 {7*l} 2* /
/

-- Set Porosities per layer
EQUALS  -- <== PARAMETER TO TUNE (Obtained using the excel file provided)
"""
            for i, porosity in enumerate(porosities, start=1):
                input_data += f"    PORO    {porosity}    4*    {odd_sequence(i)}    {even_sequence(i)}    /\n"

            input_data += """
/

-- Set Permeability values per layer
EQUALS  -- <== PARAMETER TO GET FROM POROSITY USING: PERMX = 10^(15.6*poro - 0.9)
"""
            for i, porosity in enumerate(porosities, start=1):
                input_data += f"    PERMX    {10**(15.6*porosity-0.9)}    4*    {odd_sequence(i)}    {even_sequence(i)}    /\n"

            input_data += """
/

-- Copy values to other perm directions
COPY
PERMX PERMY /
PERMX PERMZ /
/

-- Apply KvKh correction
MULTIPLY
"""
            input_data += f"     PERMZ    {permz} /  -- <== PARAMETER TO TUNE\n"

            if l == 1:
                input_data += """
/

FAULTS
 FAULTA		1	1	22	22	1	12	Y-	/
 FAULTA 	2	2	21	21	1	12	Y-	/
 FAULTA 	3	3	20	20	1	12	Y-	/
 FAULTA 	4	4	19	19	1	12	Y-	/
 FAULTA 	5	5	18	18	1	12	Y-	/
 FAULTA 	6	6	17	17	1	12	Y-	/
 FAULTA 	7	7	16	16	1	12	Y-	/
 FAULTA 	8	8	15	15	1	12	Y-	/
 FAULTA 	9	9	14	14	1	12	Y-	/
 FAULTA 	10	10	13	13	1	12	Y-	/
 FAULTA 	11	11	12	12	1	12	Y-	/
 FAULTA 	12	12	11	11	1	12	Y-	/
 FAULTA 	13	13	10	10	1	12	Y-	/
 FAULTA 	14	14	9	9	1	12	Y-	/
 FAULTA 	15	15	8	8	1	12	Y-	/
 FAULTA 	16	16	7	7	1	12	Y-	/
 FAULTA 	17	17	6	6	1	12	Y-	/
 FAULTA 	18	18	5	5	1	12	Y-	/
 FAULTA 	19	19	4	4	1	12	Y-	/
 FAULTA 	20	20	3	3	1	12	Y-	/
 FAULTA 	21	21	2	2	1	12	Y-	/
 FAULTA 	22	22	1	1	1	12	Y-	/
								
 FAULTA 	1	1	22	22	1	12	X-	/
 FAULTA 	2	2	21	21	1	12	X-	/
 FAULTA 	3	3	20	20	1	12	X-	/
 FAULTA 	4	4	19	19	1	12	X-	/
 FAULTA 	5	5	18	18	1	12	X-	/
 FAULTA 	6	6	17	17	1	12	X-	/
 FAULTA 	7	7	16	16	1	12	X-	/
 FAULTA 	8	8	15	15	1	12	X-	/
 FAULTA 	9	9	14	14	1	12	X-	/
 FAULTA 	10	10	13	13	1	12	X-	/
 FAULTA 	11	11	12	12	1	12	X-	/
 FAULTA 	12	12	11	11	1	12	X-	/
 FAULTA 	13	13	10	10	1	12	X-	/
 FAULTA 	14	14	9	9	1	12	X-	/
 FAULTA 	15	15	8	8	1	12	X-	/
 FAULTA 	16	16	7	7	1	12	X-	/
 FAULTA 	17	17	6	6	1	12	X-	/
 FAULTA 	18	18	5	5	1	12	X-	/
 FAULTA 	19	19	4	4	1	12	X-	/
 FAULTA 	20	20	3	3	1	12	X-	/
 FAULTA 	21	21	2	2	1	12	X-	/
 FAULTA 	22	22	1	1	1	12	X-	/


 FAULTB 	17	17	2	2	1	12	X-	/
 FAULTB 	18	18	4	4	1	12	X-	/
 FAULTB 	19	19	6	6	1	12	X-	/
 FAULTB 	20	20	8	8	1	12	X-	/
 FAULTB 	21	21	10	10	1	12	X-	/
 FAULTB 	22	22	12	12	1	12	X-	/
 FAULTB 	23	23	14	14	1	12	X-	/
 FAULTB 	24	24	16	16	1	12	X-	/
 FAULTB 	25	25	18	18	1	12	X-	/
 FAULTB 	26	26	20	20	1	12	X-	/
 FAULTB 	27	27	22	22	1	12	X-	/
 FAULTB 	28	28	24	24	1	12	X-	/
 FAULTB 	29	29	26	26	1	12	X-	/
 FAULTB 	30	30	28	28	1	12	X-	/
 
 
 
 FAULTB 	16	16	1	1	1	12	Y	/
 FAULTB 	17	17	3	3	1	12	Y	/
 FAULTB 	18	18	5	5	1	12	Y	/
 FAULTB 	19	19	7	7	1	12	Y	/
 FAULTB 	20	20	9	9	1	12	Y	/
 FAULTB 	21	21	11	11	1	12	Y	/
 FAULTB 	22	22	13	13	1	12	Y	/
 FAULTB 	23	23	15	15	1	12	Y	/
 FAULTB 	24	24	17	17	1	12	Y	/
 FAULTB 	25	25	19	19	1	12	Y	/
 FAULTB 	26	26	21	21	1	12	Y	/
 FAULTB 	27	27	23	23	1	12	Y	/
 FAULTB 	28	28	25	25	1	12	Y	/
 FAULTB 	29	29	27	27	1	12	Y	/
 FAULTB 	30	30	29	29	1	12	Y	/
 
 FAULTB 	16	16	1	1	1	12	X-	/
 FAULTB 	17	17	3	3	1	12	X-	/
 FAULTB 	18	18	5	5	1	12	X-	/
 FAULTB 	19	19	7	7	1	12	X-	/
 FAULTB 	20	20	9	9	1	12	X-	/
 FAULTB 	21	21	11	11	1	12	X-	/
 FAULTB 	22	22	13	13	1	12	X-	/
 FAULTB 	23	23	15	15	1	12	X-	/
 FAULTB 	24	24	17	17	1	12	X-	/
 FAULTB 	25	25	19	19	1	12	X-	/
 FAULTB 	26	26	21	21	1	12	X-	/
 FAULTB 	27	27	23	23	1	12	X-	/
 FAULTB 	28	28	25	25	1	12	X-	/
 FAULTB 	29	29	27	27	1	12	X-	/
 FAULTB 	30	30	29	29	1	12	X-	/

/

-- Control the transmissibility of the faults
MULTFLT
"""
            elif l == 2:
                input_data += """
/

FAULTS
 FAULTA		1	2	42	42	1	24	Y	/
 FAULTA 	3	4	40	40	1	24	Y	/
 FAULTA 	5	6	38	38	1	24	Y	/
 FAULTA 	7	8	36	36	1	24	Y	/
 FAULTA 	9	10	34	34	1	24	Y	/
 FAULTA 	11	12	32	32	1	24	Y	/
 FAULTA 	13	14	30	30	1	24	Y	/
 FAULTA 	15	16	28	28	1	24	Y	/
 FAULTA 	17	18	26	26	1	24	Y	/
 FAULTA 	19	20	24	24	1	24	Y	/
 FAULTA 	21	22	22	22	1	24	Y	/
 FAULTA 	23	24	20	20	1	24	Y	/
 FAULTA 	25	26	18	18	1	24	Y	/
 FAULTA 	27	28	16	16	1	24	Y	/
 FAULTA 	29	30	14	14	1	24	Y	/
 FAULTA 	31	32	12	12	1	24	Y	/
 FAULTA 	33	34	10	10	1	24	Y	/
 FAULTA 	35	36	8	8	1	24	Y	/
 FAULTA 	37	38	6	6	1	24	Y	/
 FAULTA 	39	40	4	4	1	24	Y	/
 FAULTA 	41	42	2	2	1	24	Y	/
								
 FAULTA 	42	42	1	2	1	24	X	/
 FAULTA 	40	40	3	4	1	24	X	/
 FAULTA 	38	38	5	6	1	24	X	/
 FAULTA 	36	36	7	8	1	24	X	/
 FAULTA 	34	34	9	10	1	24	X	/
 FAULTA 	32	32	11	12	1	24	X	/
 FAULTA 	30	30	13	14	1	24	X	/
 FAULTA 	28	28	15	16	1	24	X	/
 FAULTA 	26	26	17	18	1	24	X	/
 FAULTA 	24	24	19	20	1	24	X	/
 FAULTA 	22	22	21	22	1	24	X	/
 FAULTA 	20	20	23	24	1	24	X	/
 FAULTA 	18	18	25	26	1	24	X	/
 FAULTA 	16	16	27	28	1	24	X	/
 FAULTA 	14	14	29	30	1	24	X	/
 FAULTA 	12	12	31	32	1	24	X	/
 FAULTA 	10	10	33	34	1	24	X	/
 FAULTA 	8	8	35	36	1	24	X	/
 FAULTA 	6	6	37	38	1	24	X	/
 FAULTA 	4	4	39	40	1	24	X	/
 FAULTA 	2	2	41	42	1	24	X	/


 FAULTB 	59	60	59	59	1	24	Y-	/
 FAULTB 	57	58	55	55	1	24	Y-	/
 FAULTB 	55	56	51	51	1	24	Y-	/
 FAULTB 	53	54	47	47	1	24	Y-	/
 FAULTB 	51	52	43	43	1	24	Y-	/
 FAULTB 	49	50	39	39	1	24	Y-	/
 FAULTB 	47	48	35	35	1	24	Y-	/
 FAULTB 	45	46	31	31	1	24	Y-	/
 FAULTB 	43	44	27	27	1	24	Y-	/
 FAULTB 	41	42	23	23	1	24	Y-	/
 FAULTB 	39	40	19	19	1	24	Y-	/
 FAULTB 	37	38	15	15	1	24	Y-	/
 FAULTB 	35	36	11	11	1	24	Y-	/
 FAULTB 	33	34	7	7	1	24	Y-	/
 FAULTB 	31	32	3	3	1	24	Y-	/
 
 FAULTB 	58	58	55	58	1	24	X	/
 FAULTB 	56	56	51	54	1	24	X	/
 FAULTB 	54	54	47	50	1	24	X	/
 FAULTB 	52	52	43	46	1	24	X	/
 FAULTB 	50	50	39	42	1	24	X	/
 FAULTB 	48	48	35	38	1	24	X	/
 FAULTB 	46	46	31	34	1	24	X	/
 FAULTB 	44	44	27	30	1	24	X	/
 FAULTB 	42	42	23	26	1	24	X	/
 FAULTB 	40	40	19	22	1	24	X	/
 FAULTB 	38	38	15	18	1	24	X	/
 FAULTB 	36	36	11	14	1	24	X	/
 FAULTB 	34	34	7	10	1	24	X	/
 FAULTB 	32	32	3	6	1	24	X	/
 FAULTB 	30	30	1	4	1	24	X	/
 

/

-- Control the transmissibility of the faults
MULTFLT
"""

            else:
                raise Exception("Invalid level")

            input_data += f"FAULTA    {faulta} / -- Upper fault -- <== PARAMETER TO TUNE\n"
            input_data += f"FAULTB    {faultb} / -- Lower fault -- <== PARAMETER TO TUNE\n"

            dpcf = 0.0 if l == 1 else 0.02
            input_data += f"""
/

-- Dykstra-Parsons coefficient to add some heterogeneity
dpcf
{dpcf}  42 /
"""

            input_file = self.pflotran_folder + "ccs_3df.grdecl"
            with open(input_file, "w") as f:
                f.write(input_data)
            logging.info("Input file written successfully.")
        except Exception as e:
            logging.error(f"Error preparing input file: {e}")
            raise

    def supports_evaluate(self):
        return True


# Instantiate the model
pflotran_model = PFLOTRANModel()

port = int(os.getenv("PORT"))
# Serve the model using UM-Bridge
umbridge.serve_models([pflotran_model], port)
