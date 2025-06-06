import umbridge
import subprocess
import os
import logging
import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO)

class PFLOTRANModel(umbridge.Model):
    def __init__(self):
        super().__init__("pflotran_simulation")
        self.simulation_script = "/pflotran_ogs_1.8/macros/pft1.8"
        
        self.input_file_path = "/coarse_model/ccs_3df.grdecl"
        self.simulation_input = "/coarse_model/Coarse_CCS_3DF.in"
        
        self.output_folder = "/output/"
        self.output_file_path = "/coarse_model/Coarse_CCS_3DF-mas.dat"
        self.output_csv_path = "/coarse_model/Coarse_CCS_3DF.csv"
        self.cnt = 0 # Count simulations

    def get_input_sizes(self, config):
        return [12]  # 12 parameters for each layer

    def get_output_sizes(self, config):
        return [18]

    def __call__(self, parameters, config):
        try:
            # Prepare the input file for PFLOTRAN based on parameters
            self.prepare_input_file(parameters[0])
            logging.info("Input file prepared successfully.")

            # Run the PFLOTRAN simulation
            self.run_simulation()
            logging.info("Simulation ran successfully.")

            # Process the output file to extract results
            result = self.process_output_file()
            logging.info("Output processed successfully.")

            return result
        except Exception as e:
            logging.error(f"Error during simulation: {e}")
            raise

    def process_output_file(self):
        # Converts a -mas.dat into a .csv file
        replacements = {' "':",",'"':'','  ':',',' -':',-'}

        infile = self.output_file_path
        outfile= self.output_csv_path#+"_"+str(self.cnt)+".csv"

        with open(infile) as infile, open(outfile, 'w') as outfile:
            for line in infile:
                for src, target in replacements.items():
                    line = line.replace(src, target)
                outfile.write(line[1:])

        os.system("cp " +self.output_csv_path+" "+ self.output_folder+"outputfile_"+str(self.cnt))
        # Extract feature from .csv file
        df=pd.read_csv(self.output_csv_path)
        result = []
        result.append(float(df["Field fgit [m^3]"].iloc[-1])) # Total gas injected
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
        return [result]


    def run_simulation(self):
        self.cnt = self.cnt+1
        os.system("cp " +self.simulation_input+" "+ self.output_folder+"input_"+str(self.cnt))
        os.system("cp " +self.input_file_path+" "+ self.output_folder+"inputfile_"+str(self.cnt))
        try:
            os.system(self.simulation_script + " "+  self.simulation_input )

        except Exception as e:
            logging.error(f"Error running simulation: {e}")
            raise

    def prepare_input_file(self, porosities):
        try:
            input_data = f"""DIMENS
30 30 12 /

external_file ccs_3df_geom.grdecl

-- Control zone
EQUALS
FIPNOGO 5 1 7 1 7 2* /
/

-- Set Porosities per layer
EQUALS  -- <== PARAMETER TO TUNE (Obtained using the excel file provided)
"""
            for i, porosity in enumerate(porosities, start=1):
                input_data += f"    PORO    {porosity:.2f}    4*    {i}    {i}    /\n"

            input_data += """
/

-- Set Permeability values per layer
EQUALS  -- <== PARAMETER TO GET FROM POROSITY USING: PERMX = 10^(15.6*poro - 0.9)
	PERMX	212.57	4*	1	1	/
	PERMX	42.88	4*	2	2	/
	PERMX	11202.12	4*	3	3	/
	PERMX	731.09	4*	4	4	/
	PERMX	44078.84	4*	5	5	/
	PERMX	27.61	4*	6	6	/
	PERMX	0.31	4*	7	7	/
	PERMX	13.30	4*	8	8	/
	PERMX	67.10	4*	9	9	/
	PERMX	117.60	4*	10	10	/
	PERMX	4.69	4*	11	11	/
	PERMX	25.52	4*	12	12	/
/

-- Copy values to other perm directions
COPY
PERMX PERMY /
PERMX PERMZ /
/

-- Apply KvKh correction
MULTIPLY
PERMZ 0.1 /  -- <== PARAMETER TO TUNE
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
FAULTA 1.0 / -- Upper fault -- <== PARAMETER TO TUNE
FAULTB 1.0 / -- Lower fault -- <== PARAMETER TO TUNE
/

-- Dykstra-Parsons coefficient to add some heterogeneity
dpcf
0.02  42 /
"""

            with open(self.input_file_path, "w") as f:
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
