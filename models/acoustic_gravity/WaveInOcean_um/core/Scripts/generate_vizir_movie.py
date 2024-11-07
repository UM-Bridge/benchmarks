import os
import argparse
import subprocess
import re


def main(args):
    if args.solution is None:
        raise ValueError("Please provide a solution prefix name using --solution.")

    if args.result_dir is None:
        raise ValueError("Please provide a result directory using --result_dir.")

    solution_prefix = args.solution
    result_dir = args.result_dir

    if not os.path.exists(result_dir):
        raise ValueError("The result directory " + result_dir + " does not exist.")
    
    # Get Solution files and the mesh.
    mesh_file = []
    solution_files = []
    for root, _, files in os.walk(result_dir):
        for file in files:
            if file.endswith(".mesh"):
                mesh_file.append(file)
            if file.endswith(".sol") and file.startswith(solution_prefix + "."):
                solution_files.append(file)

    if not solution_files:
        raise ValueError("The solution files list for the prefix " + solution_prefix + " is empty.")

    sorted_solution_files = sorted(solution_files, key=lambda s: [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)])

    if not len(mesh_file) == 1:
        raise ValueError("There is more than one mesh file in the results directory.")
    mesh_file = mesh_file[0]

    vizir_output_folder = "VizirOut"
    if not os.path.exists(os.path.join(result_dir, vizir_output_folder)):
        os.mkdir(os.path.join(result_dir, vizir_output_folder))

    # ViZiR processing: writes binary version of the solutions files (and the mesh) while reconstructing the solution at the missing surface elements.
    sol_list_file = open(os.path.join(result_dir, solution_prefix) + ".sols", "w")
    vizir_movie_file = open(os.path.join(result_dir, "vizir.movie"), "w")
    for solution in sorted_solution_files:
        os.system("/Applications/vizir4.app/Contents/MacOS/vizir4 -checksurf -in " +  os.path.join(result_dir, mesh_file) +  " -sol " + os.path.join(result_dir, solution))
        # subprocess.run(["/Applications/vizir4.app/Contents/MacOS/vizir4", " -checksurf -in ", os.path.join(result_dir, mesh_file), " -sol ", os.path.join(result_dir, solution)], stdout=subprocess.PIPE)
        solution_insert_index = solution.find(".sol")
        mesh_insert_index = mesh_file.find(".mesh")
        binary_solution_name = solution[:solution_insert_index] + ".surfOK" + solution[solution_insert_index:] + "b"
        binary_mesh_name = mesh_file[:mesh_insert_index] + ".surfOK" + mesh_file[mesh_insert_index:] + "b"
        sol_list_file.write(binary_solution_name + "\n")
        vizir_movie_file.write(binary_mesh_name + " " + binary_solution_name + " " + vizir_output_folder + "/" + solution[:solution_insert_index] + ".jpg\n")

    sol_list_file.close()
    vizir_movie_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
      "--solution",
      "-f",
      "-sol",
      help = "prefix for the solution file",
      type = str,
      )
    
    parser.add_argument(
      "--result_dir",
      help = "Path to the result dir",
      type = str,
      )
  
    args = parser.parse_args()
    main(args)
    