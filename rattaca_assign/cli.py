#!/usr/bin/env python3
'''
Command-line interface for RATTACA assignment tool.

Usage: rattaca_assign -c colony_dataframe.csv -p predictions.csv -o /output/path/output_prefix -r request1.json request2.json...
'''

import argparse
import sys
from rattaca_assign.core.assign import run_assignments
from rattaca_assign.core.model_utils import which_requests


def parse_args(args=None):

    '''Parse command line arguments for the RATTACA tool.'''
    
    if args is None:
        args = sys.argv[1:]
        
    parser = argparse.ArgumentParser(
        prog='RATTACA assignment',
        description='''Uses permutation to maximize per-project
            between-group differences in predicted traits when
            assigning rats to new projects''')
            
    parser.add_argument('-o', '--output_prefix', 
        required=False, type=str, nargs=1, 
        dest='output_prefix', help='<output file prefix: /desired/path/plus_file_prefix>')
    parser.add_argument('-c', '--colony_df', 
        required=True, type=str, nargs=1, 
        dest='colony_dataframe', help='<HS West colony dataframe csv>')
    parser.add_argument('-p', '--predictions', 
        required=False, type=str, nargs=1, 
        dest='predictions', help='<rattaca predictions csv')
    parser.add_argument('-r', '--request_files', 
        required=True, type=str, 
        nargs='+', dest='requests', help='<requests jsons')
    parser.add_argument('-e', '--exclude_rfids', 
        required=False, type=str, nargs=1,
        dest='exclude', help='RFIDs in colony data to exclude from assignment')

    return parser.parse_args(args)
