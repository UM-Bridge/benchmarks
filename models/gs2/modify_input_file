#!/usr/bin/env python3

''' This script modifies selected parameters of a MOOSE/Achlys input file
according to a separate JSON "modifier" file. The JSON file must have the 
same hierachical structure as the target MOOSE input file.
'''

__author__ = 'Mikkel Bue Lykkegaard'
__copyright__ = 'Copyright 2022, digiLab'
__credits__ = ['Mikkel Bue Lykkegaard', 'Tim Dodwell', 'Freddy Wordingham']
__license__ = 'MIT'
__version__ = '0.1'
__maintainer__ = 'Mikkel Bue Lykkegaard'
__email__ = 'mikkel@digilab.co.uk'
__status__ = "Development"

# import core packages
import collections.abc as collections

# import MOOSE packages
import pyhit
import moosetree

def flatten_dict(d, parent_key='', sep='/'):
    
    '''Helper function that turns a nested dict into a flattened dict.
    Adapted from https://stackoverflow.com/questions/6027558.
    
    Parameters
    ----------
    d : dict
        Nested input dict.
    parent_key : str, optional
        Additional key to append in front of every key (used in recursion).
    sep : str, optional
        String to use as separator of the flattened keys. Default is '/'.
    '''
    
    # set up empty list to hold the dict items.
    items = []
    
    # iterate through the items of the current level.
    for k, v in d.items():
        
        # construct a string for the flattened key.
        #new_key = parent_key + sep + k if parent_key else k
        # NB: This line has been modified from the original to add an 
        # additional separator at the beginning for the (flattened) key.
        new_key = parent_key + sep + k if parent_key else sep + k
        
        # if the value of an item is another dict, call this function 
        # recursively and extend.
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
            
        # otherwise, append the key and value pair.
        else:
            items.append((new_key, v))
            
    # turn the list into a dict and return.
    return dict(items)

def modify_parameter(root, parameter_path, new_value):
    
    '''Modifies a parameter in the MOOSE input file
    
    Parameters
    ----------
    root : pyhit.Node
        The root node of the MOOSE input file.
    parameter_path : str
        The full path of the parameter to be modified.
    new_value : float
        The new parameter value.
    '''
    
    # split the parameter path into the path of the parameter node
    # and the parameter name.
    path, _, parameter = parameter_path.rpartition('/')
    
    # find the node in the MOOSE tree.
    node = moosetree.find(root, func=lambda n: n.fullpath == path)
    
    # set the value.
    node[parameter] = new_value

# do this stuff if called from the command line.
if __name__ == '__main__':
    
    # import core packages.
    import argparse
    import json

    # set up a parser.
    parser = argparse.ArgumentParser(description='Update MOOSE input file.')
    
    # add the modifier file to the parser.
    parser.add_argument('modifier', metavar='Modifier', type=str,
                        help='A JSON file containing the modified parameter names and values')
    # add the input file.
    parser.add_argument('infile', metavar='Input file', type=str,
                        help='A MOOSE input file to be modified')
    # add the (optional) output file.
    parser.add_argument('-o', '--outfile', dest='outfile', type=str, default=None,
                        help='An output file name (default: same as input file (overwrite))')
    
    # parse the args from command line.
    args = parser.parse_args()
    
    # set the outfile to infile if no outfile was supplied.
    if args.outfile is None:
        args.outfile = args.infile
    
    # load the modifier JSON and flatten the dict.
    with open(args.modifier) as modifier_file:
        modifier = flatten_dict(json.load(modifier_file)['Achlys'])
    
    # read the MOOSE input file
    root = pyhit.load(args.infile)
    
    # modify each parameter in the modifier dict.
    for k, v in modifier.items():
        modify_parameter(root, k, v)
    
    # Write the modified file
    pyhit.write(args.outfile, root)
