import subprocess
import os
import shutil
import sys

from substitute_env_variable import SubstituteEnvVariable


class ComputeCodeStats:
    """
    Use the `cloc` (https://github.com/AlDanial/cloc/releases) program several times to generate statistics
    for each module of OndoMathX.
    
    Consolidated stats are generated at the end, by amalgamating the module reports. There are two versions: one with
    everything and one wihout the tests and model instances.
    """
    
    def __init__(self, output_directory):
        
        '''
        Compute the statistics about OndoMathX.
        \param output_directory Directory into which output and work files are written.
        '''

        self.__CheckCloc()

        self.__output_directory = SubstituteEnvVariable(output_directory)
    
        if not os.path.exists(self.__output_directory):
            os.makedirs(self.__output_directory)
    
        # OndoMathX root folder is twice removed from the one this script is found.
        self.__ondomathx_root_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
        self.__ondomathx_root_directory = os.path.normpath(self.__ondomathx_root_directory)

        folder_to_amalgamate = []

        folder_to_amalgamate.append(self.__ComputeStatsForSubdir(".", do_iterate_through_subdirs = False, output_name = "Root"))
        folder_to_amalgamate.append(self.__ComputeStatsForSubdir("Scripts"))
        folder_to_amalgamate.append(self.__ComputeStatsForSubdir("Sources", generate_markdown=True))

        # Write amalgamations of reports
        self.__Amalgamate(folder_to_amalgamate, "reports_all")
        

    def __CheckCloc(self):
        """Check cloc is installed on the computer, and if not suggest to install it"""
        script_name = "cloc"

        ret = shutil.which(script_name)

        if not ret:
            raise Exception("'{}' is not installed on your system; please install it so that the executable is in your path! Several installers are available, please check https://github.com/AlDanial/cloc/blob/master/cloc.".format(script_name))
        
        return ret


    def __ComputeStatsForSubdir(self, subdirectory, do_iterate_through_subdirs = True, output_name = None, exclude_files = None, generate_markdown = False):
        """
        Compute the stats for the sources in 'subdirectory'; they will be written in a 'output_name' subdirectory of the output directory.

        \param[in] do_iterate_through_subdirs If true, will recursively iterate through all subdirectories of 'subdirectory'
        \param[in] exclude_files The paths of the files to ignore. Should be added relatively to 'subdirectory'.
        \param[in] output_name If defined, results will be written in self.__output_directory / 'output_name'.
        If not, self.__output_directory / 'subdirectory' will be used instead
        \param[in] generate_markdown If True, also generate Markdown report 

        \return Output directory into which results were written.
        """
        print("======= COMPUTE STATISTICS FOR {} =======".format(subdirectory))
        
        if output_name:
            output_directory = os.path.join(self.__output_directory, output_name)
        else:
            output_directory = os.path.join(self.__output_directory, subdirectory)
        
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        
        cmd = ["cloc", \
               os.path.join(self.__ondomathx_root_directory, subdirectory), \
               "--force-lang=C++,hxx", \
               "--force-lang=C++,doxygen", \
               "--force-lang=CMake,cmake", \
               "--exclude-ext=scl,geo,sol,lua", \
                "--fullpath ", \
                "--not-match-d='htmlcov'"
        ]

        if not do_iterate_through_subdirs:
            cmd.append("--no-recurse")

        if exclude_files:            
            exclude_files_path = os.path.join(output_directory, "input_exclude_list.txt")

            with open(exclude_files_path, "w") as excluded_file_stream:
                for elem in exclude_files:

                    to_ignore_file = os.path.normpath(os.path.join(self.__ondomathx_root_directory, subdirectory, elem))
                    excluded_file_stream.write("{}\n".format(to_ignore_file))

            cmd.extend(('--exclude-list-file', exclude_files_path))
         
        cmd.extend((
                '--ignored', \
                os.path.join(output_directory, 'ignored.txt'), \
                '--counted', \
                os.path.join(output_directory, 'files_considered.txt'), \
                "--report-file ", \
                os.path.join(output_directory, 'report.txt')))

        subprocess.Popen(' '.join(cmd), shell = True).communicate()

        if generate_markdown:
            del cmd[-1]
            cmd.extend((os.path.join(output_directory, 'report.md'), "--md"))

            subprocess.Popen(' '.join(cmd), shell = True).communicate()
        

        return output_directory
        
    def __Amalgamate(self, list_dir, report_name):
        """Amalgamate several report files into one.
        
        \param[in] list_dir List of--exclude-lang= subdirectories of self.__ondomathx_src for which we track the reports. 
        Only the subdirectory name is expected, not the full path (e.g. 'Utilities').
        \param[in] report_name Name of the report file. It will be written in self.__output_directory.
        """
        report_list = []
        for subdir in list_dir:
            report_list.append(os.path.join(self.__output_directory, subdir, "report.txt"))
        
        cmd = ("cloc", \
               "--sum-reports ", \
               " ".join(report_list), \
               "--report-file ", \
               os.path.join(self.__output_directory, report_name + ".txt")
               )
            
        subprocess.Popen(' '.join(cmd), shell = True).communicate()

        # Do it again in Markdown format (easier to transform into a pdf)
        cmd = ("cloc", \
               "--sum-reports ", \
               " ".join(report_list), \
               "--md", \
               "--report-file ", \
               os.path.join(self.__output_directory, report_name + ".md")
               )
            
        subprocess.Popen(' '.join(cmd), shell = True).communicate()



if __name__ == "__main__":    
    ComputeCodeStats('./CodeStats')
