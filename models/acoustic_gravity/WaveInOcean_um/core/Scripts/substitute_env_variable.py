import os


def SubstituteEnvVariable(directory, begin_separator = "${", end_separator = "}"):
    """Replace environment variables by their value (which must be properly defined).

    Syntax to provide is:
        ${ENVIRONMENT_VARIABLE}
    where {ENVIRONMENT_VARIABLE} is the name of the environment variable (e.g. ${HOME}).

    Syntax may actually be customized with begin_separator and end_separator parameters.

    There is an assumption here it is a path we are figuring out (os.path.normpath is used).

    \return Name of the directory with the value of the environment variable substituted.
    """
    
    begin_index = directory.find(begin_separator)
    
    while begin_index != -1:
                
        try:
            end_index = directory.index(end_separator)
        except ValueError:
            raise Exception("There are not enough end separator in SubstituteEnvVariable() (typically "
            "you have a \{ for which there are no matching \} with default separators).")
            
        name_env_variable = directory[begin_index + len(begin_separator) : end_index]
        env_variable = directory[begin_index : end_index + len(end_separator)]
        
        directory = directory.replace(env_variable, os.environ[name_env_variable])
            
        begin_index = directory.find(begin_separator)

    return directory
    
    
    
if __name__ == "__main__":
    
    os.environ["FOO"] = "foo"
    os.environ["BAR"] = "bar"
    os.environ["BAZ"] = "baz"
    
    directory = "/truc/bidule/${FOO}/subdir_${BAR}/${BAZ}/end"
    expected = "/truc/bidule/foo/subdir_bar/baz/end"
    
    SubstituteEnvVariable(directory)
    
    if SubstituteEnvVariable(directory) != expected:
       raise Exception("|{}| was expected but we got \n|{}|".format(expected, SubstituteEnvVariable(directory)))
    
    directory = "/truc/bidule/${===FOO===}/subdir_${===BAR===}/${===BAZ===}/end"
    
    more_complex = SubstituteEnvVariable(directory, begin_separator = "${===", end_separator = "===}") 
    if more_complex != expected:
        raise Exception("|{}| was expected but we got \n|{}|".format(expected, more_complex))
    