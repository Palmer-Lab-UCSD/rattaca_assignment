#!/usr/bin/env python3
'''
Command-line interface for RATTACA assignment tool.

Usage: rattaca_assign -c colony_dataframe.csv -p predictions.csv -o /output/path/output_prefix -r request1.json request2.json...
'''

import argparse
import sys
# from rattaca_assign.core.assign import run_assignments
# from rattaca_assign.core.model_utils import which_requests


def parse_args(args=None):

    '''Parse command line arguments for the RATTACA tool.'''
    
    if args is None:
        args = sys.argv[1:]
        
    parser = argparse.ArgumentParser(
        prog='RATTACA assignment',
        description='''Uses permutation to maximize per-project
            between-group differences in predicted traits when
            assigning rats to new projects''')

    # positional arguments
    parser.add_argument('step', 
        choices=['new_assignments','update_assignments','rattaca_results'],
        type=str, 
        help='the step to execute: [new_assignments, update_assignments, rattaca_results]')
        
    # non-positional arguments
    parser.add_argument('-o', '--output_dir', 
        required=False, type=str, nargs=1, 
        dest='output_dir', help='<output directory: /desired/output/folder>')
    parser.add_argument('-f', '--output_prefix', 
        required=False, type=str, nargs=1, 
        dest='output_prefix', help='<output file prefix>')
    parser.add_argument('-c', '--colony_df', 
        required=False, type=str, nargs=1, # originally True
        dest='colony_dataframe', help='<HS West colony dataframe csv>')
    parser.add_argument('-p', '--predictions', 
        required=False, type=str, nargs=1, 
        dest='predictions', help='<rattaca predictions csv>')
    parser.add_argument('-s', '--predictions_summary', 
        required=False, type=str, nargs=1, 
        dest='preds_summary', help='<rattaca predictions summary csv>')
    parser.add_argument('-d', '--preds_dir', 
        required=False, type=str, nargs=1, 
        dest='preds_dir', help='path housing all rattaca predictions')
    parser.add_argument('-r', '--request_files', 
        required=False, type=str, # originally True
        nargs='+', dest='requests', help='<requests jsons>')
    parser.add_argument('-e', '--exclude_rfids', 
        required=False, type=str, nargs=1,
        dest='exclude', help='RFIDs in colony data to exclude from assignment')
    parser.add_argument('-m', '--request_map', 
        required=False, type=str, nargs=1,
        dest='request_map', help='path to request map csv (rows=request, cols=traits)') 
    parser.add_argument('-a', '--proposed_assignments',
        required=False, type=str, nargs=1,
        dest='assignments', help='path to csv of previously proposed assignments')
    parser.add_argument('-u', '--update_requests',
        required=False, type=str, nargs=1,
        dest='updates', help='path to json file with requests to be updated by shipping sheets')
    parser.add_argument('-v', '--verbose',
            action='store_true',
            dest='verbose', help='enable verbose output for debugging')
    
    return parser.parse_args(args)
