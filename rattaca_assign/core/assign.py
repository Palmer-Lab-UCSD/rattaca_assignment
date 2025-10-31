'''
This is RATTACA's core assignment algorithm.

This module contains the main assignment logic that distributes rats to 
different projects based on their genetic predictions and other criteria, plus
utility functions to format and save assignment results.
'''

import sys
import pandas as pd
import json
import glob
from random import choice
from pprint import pprint
from os import path
from pathlib import Path
from datetime import datetime
from itertools import permutations
from rattaca_assign.core.model_utils import * 
from rattaca_assign.core.base_utils import *

date_today = datetime.today().date() 
datestamp = date_today.strftime('%Y%m%d')

def run_assignments(args):
    '''
    Run the assignment algorithm to distribute rats to projects.

    This is the top-level assignment function that integrates request-specific 
    assignment methods to distribute rats among all projects included in args.
    
    Args:
        args: command line arguments containing request files.
        
    Returns:
        A dictionary mapping projects to assigned RFIDs.
    '''
    # get request types
    requests_made = which_requests(args)
    
    # identify which types of requests have been made
    breeders_requested = 'hsw_breeders' in requests_made
    rattaca_requested = 'rattaca' in requests_made
    random_requested = 'random' in requests_made

    # load and group request objects
    request_objects = load_request_files(args)
    
    # Get request objects by type
    rattaca_requests = request_objects['rattaca'] if rattaca_requested else []
    random_requests = request_objects['random'] if random_requested else []
    breeder_request = request_objects['hsw_breeders'][0] if breeders_requested and request_objects.get('hsw_breeders') else None
    
    # Initialize collections
    all_requests = []
    non_breeder_requests = []
    non_rattaca_requests = []
    
    # Add breeder request if it exists
    if breeders_requested and breeder_request:
        all_requests.append(breeder_request)
        
    # Add RATTACA requests
    if rattaca_requested and rattaca_requests:
        all_requests.extend(rattaca_requests)
        non_breeder_requests.extend(rattaca_requests)
        
    # Add random requests
    if random_requested and random_requests:
        all_requests.extend(random_requests)
        non_breeder_requests.extend(random_requests)
        non_rattaca_requests.extend(random_requests)

    # STAGE 1: priority assignment to HSW breeders
    if breeders_requested and breeder_request:
        priority_breeders = prioritize_breeders(
            breeder_request=breeder_request,
            non_breeder_requests=non_breeder_requests
        )
        print(f'priority_breeders: {priority_breeders}')

        # assign breeders by priority
        breeder_request.assign_hsw_breeders(
            priority_breeders,
            non_breeder_requests=non_breeder_requests)
        
        # remove assigned breeders from availability for other projects
        for request in non_breeder_requests:
            request.remove(priority_breeders)
    
    # STAGE 2: RATTACA assignments
    if rattaca_requested and rattaca_requests:
        # continue assigning RATTACA rats until all RATTACA requests are satisfied
        rattaca_assignment(
            rattaca_requests=rattaca_requests, 
            non_rattaca_requests=non_rattaca_requests)  # Fixed the syntax error here
    
    # STAGE 3: Complete breeder assignments after RATTACA
    if breeders_requested and breeder_request:
        breeder_requests = [breeder_request]  # Create a list for consistency with the rest of the code
        if not breeder_request.is_satisfied():
            # Find any new priority breeders after RATTACA assignments
            new_priority_breeders = prioritize_breeders(
                breeder_request=breeder_request,
                non_breeder_requests=non_breeder_requests
            )
            if new_priority_breeders:
                breeder_request.assign_manual_hsw_breeders(new_priority_breeders)
                # Remove from random projects
                for request in random_requests:
                    request.remove(new_priority_breeders)
            
            # For families that have contributed to RATTACA, ensure a male and female 
            # from each family goes to breeding if available
            # complete_breeding_assignments(breeder_request) - defunct function, fix this section
    
    # STAGE 4: Only now assign to random projects
    if random_requested and random_requests:
        # At this point, all RATTACA and priority breeder assignments are complete
        # Random projects get what's left
        for random_req in random_requests:
            # Ensure random projects only use rats from families that have already
            # contributed to RATTACA or breeders when possible
            assign_to_random_projects(random_req, rattaca_requests, breeder_requests)
    
    # return all request objects for output processing
    result_requests = []
    if rattaca_requested and rattaca_requests:
        result_requests.extend(rattaca_requests)
    if breeders_requested and breeder_request:
        result_requests.append(breeder_request)
    if random_requested and random_requests:
        result_requests.extend(random_requests)
    
    return result_requests


def rattaca_assignment(rattaca_requests, non_rattaca_requests = None):
    '''
    Run RATTACA assignment algorithm until all requests are satisfied.
    
    Args:
        rattaca_requests: List of RATTACA request objects
        non_rattaca_requests: List of non-RATTACA request objects. Availability
        for these requests will be automatically updated as RFIDs are assigned 
        to RATTACA requests.
    '''

    assignment_round = 0
    while any(not request.is_satisfied() for request in rattaca_requests):
        assignment_round += 1
        print(f'RATTACA ASSIGNMENT ROUND {assignment_round}')
        
        open_requests = [req for req in rattaca_requests if not req.is_satisfied()]
        if not open_requests:
            break
            
        # run one round of RATTACA permutation assginment
        assigned_rfids = permute_one_round(open_requests)
        
        # after each round, update breeder requests to reflect the new assignments
        if non_rattaca_requests and assigned_rfids:
            for request in non_rattaca_requests:
                request.remove(assigned_rfids)

def permute_rattaca(all_requests):
    '''
    Assign rats to RATTACA projects using permutation to maximize overall delta.
    
    Args:
        all_requests: a list of RATTACA request objects to process.
        
    Returns:
        A dictionary mapping projects to assigned RFIDs.
    '''
    assignment_results = {}
    
    # continue until all requests are satisfied
    assignment_round = 0
    while any(not request.is_satisfied() for request in all_requests):
        # get unsatisfied requests
        open_requests = \
            [request for request in all_requests if not request.is_satisfied()]
        n_open_requests = len(open_requests)
        
        # initialize variables to store stats from the eventual best permutation
        best_delta = 0
        best_permutation = None
        best_rfids = None

        assignment_round += 1
        print(f'ASSIGNMENT ROUND {assignment_round}')
        
        # try all possible permutations of request orders
        for permuted_order in permutations(range(n_open_requests)):
            proposed_rfids = {}
            permutation_delta = 0
            project_deltas = []

            # for each project request in the current permutation... 
            for current_request_idx in permuted_order:
                current_request = open_requests[current_request_idx]
                # ...propose animals for assignment to the project and add
                # the project's delta to the total delta for the permutation
                project_delta, proposal_rfids = \
                    current_request.proposal(proposed_rfids) 
                permutation_delta += project_delta
                project_deltas.append(project_delta)
                proposed_rfids[current_request_idx] = proposal_rfids

            # keep the proposed rfids for assignment from the permutation 
            # that maximizes delta
            if permutation_delta > best_delta:
                best_delta = permutation_delta
                best_permutation = permuted_order
                best_rfids = proposed_rfids

        # assign rats to each request in the permuted order that maximizes delta
        for project_index in best_permutation:
            project_to_assign = open_requests[project_index]
            other_projects = [proj for i, proj in enumerate(open_requests) 
                              if i != project_index]
            
            assign_rfids = best_rfids[project_index]
            remove_rfids = [rfid for remaining_rfids in 
                            [rfid_val for idx_key, rfid_val in 
                             best_rfids.items() if idx_key != project_index] 
                                for rfid in remaining_rfids]

            if not project_to_assign.is_satisfied():
                project_to_assign.assign(assign_rfids)
                for rm_from_project in other_projects:
                    rm_from_project.remove(assign_rfids)
                
                # store results for this project
                project_name = project_to_assign.project
                if project_name not in assignment_results:
                    assignment_results[project_name] = []
                assignment_results[project_name].extend(assign_rfids)
    
    return all_requests


def permute_one_round(open_requests):
    '''
    Perform one round of permutation-based assignment for non-breeder requests.
    
    Args:
        open_requests: a list of unsatisfied non-breeder request objects to process.
        
    Returns:
        A list of RFIDs assigned in this round.
    '''
    n_open_requests = len(open_requests)
    if n_open_requests == 0:
        return []
    
    # Initialize variables to store stats from the eventual best permutation
    best_delta = 0
    best_permutation = None
    best_rfids = None
    
    # Try all possible permutations of request orders
    for permuted_order in permutations(range(n_open_requests)):
        proposed_rfids = {}
        permutation_delta = 0
        project_deltas = []

        # For each project request in the current permutation... 
        for current_request_idx in permuted_order:
            current_request = open_requests[current_request_idx]
            # ...propose animals for assignment to the project
            project_delta, proposal_rfids = current_request.proposal(proposed_rfids) 
            permutation_delta += project_delta
            project_deltas.append(project_delta)
            proposed_rfids[current_request_idx] = proposal_rfids

        # Keep the proposed rfids from the permutation that maximizes delta
        if permutation_delta > best_delta:
            best_delta = permutation_delta
            best_permutation = permuted_order
            best_rfids = proposed_rfids

    # Assign rats to each request in the permuted order that maximizes delta
    assigned_this_round = []
    for project_index in best_permutation:
        project_to_assign = open_requests[project_index]
        other_projects = [proj for i, proj in enumerate(open_requests) 
                          if i != project_index]
        
        assign_rfids = best_rfids[project_index]
        if assign_rfids and not project_to_assign.is_satisfied():
            project_to_assign.assign(assign_rfids)
            assigned_this_round.extend(assign_rfids)
            for rm_from_project in other_projects:
                rm_from_project.remove(assign_rfids)
    
    return assigned_this_round


def output_assignment_results(args, assigned_requests, output_prefix=None):
    '''
    Output assignment results to CSV files.
    
    Args:
        assigned_requests: List of request objects with assigned rats
        output_prefix: path and file basename prefix to write output files.
        args: Command line arguments containing colony dataframe path
        
    Returns:
        A dataframe containing all assignment results.
    '''

    # read in the colony dataframe
    df = pd.read_csv(args.colony_dataframe[0])

    # format dates
    df['dob'] = pd.to_datetime(df['dob'], format='mixed').dt.date
    df['dow'] = pd.to_datetime(df['dow'], format='mixed').dt.date
    df['dorfid'] = pd.to_datetime(df['dorfid'], format='mixed').dt.date

    # add assignment columns                       
    df['project_name'] = 'not_assigned'
    df['request_name'] = 'not_assigned'
    df['assignment'] = 'not_assigned'

    # organize columns 
    assign_cols = ['generation','animalid','accessid','rfid','earpunch','sex',
        'coatcolor','rack_number','rack_position','plate_id','well_position',
        'dob','dow','dorfid','breederpair','dam','sire','littersize','num_weaned',
        'num_sacked','project_name','request_name','assignment', 'comments']
    col_dtypes = {'generation': 'Int64', 'animalid':'str','accessid':'int', 
        'rfid':'str', 'earpunch':'str', 'sex':'str', 'coatcolor':'str', 
        'rack_number':'int','rack_position':'str','plate_id':'str', 
        'well_position':'str', 'dob':'str','dow':'str', 'dorfid':'str',
        'breederpair': 'str', 'dam':'str', 'sire':'str','littersize':'Int64', 
        'num_weaned':'Int64', 'num_sacked':'Int64', 'project_name': 'str',
        'request_name':'str', 'assignment':'str'}
    df = df[assign_cols]
    df = df.astype(col_dtypes).sort_values('animalid')
    assign_cols.remove('comments')
    
    # get assigned RFIDs from each request
    assign_names = []
    for req in assigned_requests:

        project_name = req.project
        request_name = req.request_name
        assign_name = req.assignment_name
        assign_names.append(assign_name)

        if req.assignment_type == 'rattaca':
            assigned_rfids = list(req.assigned_rats['assigned_total'].keys())
        else:
            assigned_rfids = list(req.assigned_rats.keys())
                
        # write assignments to the dataframe
        df.loc[df['rfid'].isin(assigned_rfids), 'project_name'] = project_name
        df.loc[df['rfid'].isin(assigned_rfids), 'request_name'] = request_name
        df.loc[df['rfid'].isin(assigned_rfids), 'assignment'] = assign_name
        df[assign_name] = df['rfid'].isin(assigned_rfids).astype(int)

        # save request-specific assignments to file
        req_df = df[df['assignment']==assign_name]
        req_cols = assign_cols + ['comments']
        req_df = req_df[req_cols]
        req_df.sort_values('animalid', inplace=True)
        req_file = f'{output_prefix}_assign_{assign_name}_{datestamp}.csv'
        req_df.to_csv(req_file, index=False, na_rep='')
        print(f'\n{assign_name} assignments saved to {req_file}')

    # set the assignment name order to put breeders first
    if 'hsw_breeders' in assign_names:
        assign_names.remove('hsw_breeders')
        assign_names.insert(0, 'hsw_breeders')
        
    # write all assignments to file
    assign_cols = assign_cols + assign_names + ['comments']
    df = df[assign_cols]
    assignments_file = f'{output_prefix}_assignments_{datestamp}.csv'
    df.to_csv(assignments_file, index=False, na_rep='')    
    print(f'\nAll assignments saved to {assignments_file} \n')

    return(df)


def output_assignment_preds(assignments, preds, outdir, requests, request_map=None, preds_dir=None, summary=None):
    '''
    Merge assignment results with predictions and save to CSV files formatted 
    for requesters.
    
    Args:
    
        assignments: A csv path or Pandas dataframe with assignment mapped to 
            requests, as output by output_assignment_results() or update_assignments().
    
        preds: A csv path or Pandas dataframe with trait predictions, as output
            by the 'rattaca' R package. 
        
        outdir: The directory in which to save output files.
        
        requests: A list of request file paths or request dictionaries. 
        
        request_map: Optional. A csv path with rows = requests and columns = traits.
            Used to define which trait predictions should be saved together for a 
            given request.
        
        preds_dir: Optional. The directory path holding all prediction results.
            Only used if preds is not a path string.

        summary: Optional. A csv path to a predictions summary file, as output by
            rattaca::summarize_preds().
        
    Returns:
        A dataframe containing all assignment results.
    '''
    
    if isinstance(assignments, str):
        assign_df = pd.read_csv(assignments, dtype={'rfid': str})
    else:
        assign_df = assignments
    
    if isinstance(preds, str):
        preds_df = pd.read_csv(preds, dtype={'rfid': str})
        preds_dir = os.path.dirname(preds)
    else:
        preds_df = preds
        if preds_dir is None:
            raise ValueError('Must input a preds_dir path')

    if isinstance(request_map, str):
        request_map = pd.read_csv(request_map)
    
    gen = assign_df['generation'][0]

    # format dates
    assign_df['dob'] = pd.to_datetime(assign_df['dob'], format='mixed').dt.date
    assign_df['dow'] = pd.to_datetime(assign_df['dow'], format='mixed').dt.date

    # save the generation number
    gen = assign_df['generation'].tolist()[0]
    
    # organize output columns 
    assign_cols = ['rfid','animalid','earpunch','sex','coatcolor','generation',
                   'dob','dow','project_name','request_name','assignment', 'comments']
    col_dtypes = {'generation': 'Int64', 'animalid':'str', 'rfid':'str', 
                  'earpunch':'str', 'sex':'str', 'coatcolor':'str', 'dob':'str', 'dow':'str', 
                  'project_name':'str', 'request_name':'str', 'assignment':'str','comments':'str'}

    # get request traits - handle both args and direct request inputs
    request_traits = {}
    assign_traits = {}
    assign_names = {}
    
    # read in requests from a list
    for req in requests:
    
        # read from file paths
        if isinstance(req, str):
            with open(req, 'r') as rf:
                req_data = json.load(rf)
                req_name = req_data['request_name']
                req_trait = req_data['trait']
                if req_trait is not None:
                    request_traits[req_name] = req_trait
                assign_traits[req_name] = req_data['trait']
                assign_names[req_name] = req_data['assignment_name']
                
        # read from dictionaries
        elif isinstance(req, dict):
            req_name = req['request_name']
            req_trait = req['trait']
            if req_trait is not None:
                request_traits[req_name] = req_trait
            assign_traits[req_name] = req_data['trait']
            assign_names[req_name] = req_data['assignment_name']
            
    # process each request
    for req_name in request_traits.keys():

        print(f'Processing request: {req_name} \n')
        req_trait = request_traits[req_name]

        # subset to assigned rats
        req_assignments = assign_df[assign_df['request_name'] == req_name]
        req_assignments = req_assignments.sort_values('rfid')
        req_assignments = req_assignments[assign_cols].copy()
        assign_cols_out = assign_cols.copy()
        assign_cols_out.remove('comments')
                
        # extract predictions for all requested traits (if provided)
        if request_map is not None:
            map_row = request_map[request_map['request_name'] == req_name]
            req_traits = [col for col in map_row.columns if map_row[col].iloc[0] == 1]

            # get all columns for requested traits
            all_cols = ['rfid']
            for trait in req_traits:
                all_cols.extend([trait, f'{trait}_rank', f'{trait}_zscore'])
            
            req_preds = preds_df[all_cols].copy()
            
            for req_trait in req_traits:
                group_col = f'{req_trait}_group'
                # create temporary dataframe for trait_groups function
                trait_cols = ['rfid', req_trait, f'{req_trait}_rank', f'{req_trait}_zscore']
                temp_df = req_preds[trait_cols].copy()
                
                # get assignment groups for each trait
                req_preds[group_col] = trait_groups(preds=temp_df, trait=req_trait)

            # assemble all trait results
            pred_cols = ['rfid']
            for trait in req_traits:
                pred_cols.extend([trait, f'{trait}_rank', f'{trait}_zscore', f'{trait}_group'])
            
            req_preds = req_preds[pred_cols].copy()
            pred_cols.remove('rfid')
            
        # or, extract predictions only for the assignment trait
        else:
            pred_cols = ['rfid', req_trait, f'{req_trait}_rank', f'{req_trait}_zscore']
            req_preds = preds_df[pred_cols].copy()
            group_col = f'{req_trait}_group'
            req_preds[group_col] = trait_groups(preds=req_preds, trait=req_trait, n_groups=2)
            pred_cols = pred_cols + [group_col]
            pred_cols.remove('rfid')

        # merge asssignments and predictions
        req_preds = req_assignments.merge(req_preds, on='rfid')
        out_cols = assign_cols_out + pred_cols + ['comments']
        req_preds = req_preds[out_cols]

        if Path(outdir).parts[-1] != 'request_results':
            req_outdir = os.path.join(outdir, 'request_results')
            os.makedirs(outdir, exist_ok = True)
        else:
            req_outdir = outdir
            
        req_outdir = os.path.join(req_outdir, req_name)
        os.makedirs(req_outdir, exist_ok = True)

        assign_outfile = os.path.join(req_outdir, f'{req_name}_results_BLIND.csv')
        preds_outfile = os.path.join(req_outdir, f'{req_name}_results.csv')
        
        req_assignments.to_csv(assign_outfile, index=False)
        req_preds.to_csv(preds_outfile, index=False)
        print(f'{req_name} assignments + predictions saved to {preds_outfile} \n')

        # copy accessory prediction results into assignment results
        if request_map is not None:
            traits_to_copy = req_traits
        else:
            traits_to_copy = [req_trait]

        print(f'Copying accessory predictions data for {req_name} \n')

        if summary is not None:
            if os.path.exists(summary):
                req_summary = pd.read_csv(summary)
                req_summary = req_summary[req_summary['trait'].isin(traits_to_copy)]
                sumfile = os.path.join(req_outdir, f'{req_name}_summary.csv')
                req_summary.to_csv(sumfile, index=False)
                print(f'Copied {req_name} traits from summary file {summary} \n')
            else:
                print(f'Summary file not found: {summary} \n')

        for trait in traits_to_copy:
            trait_preds_dir = os.path.join(preds_dir, trait)
            power_analyses = glob.glob(os.path.join(trait_preds_dir, f'{trait}_test_power*_summary.csv'))
            corr_dat = glob.glob(os.path.join(trait_preds_dir, f'{trait}_corr_*_traits.csv'))
            corrplot = glob.glob(os.path.join(trait_preds_dir, f'{trait}_corrplots.png'))

            print(f'Processing trait: {trait} \n')
            print(f'Power analyses found: {power_analyses}')
            print(f'Trait correlations found: {corr_dat}')
            print(f'Correlation plots found: {corrplot}')
            trait_preds_dat = power_analyses + corr_dat + corrplot

            for file in trait_preds_dat:
                if os.path.exists(file):  
                    shutil.copy(file, req_outdir)
                    print(f'Copied: {os.path.basename(file)} \n')
                else:
                    print(f'File not found: {file}')
            print(f'\n')

        # save assignment plots
        print(f'Plotting {assign_traits[req_name]} assignments for {req_name} \n')

        plot_assignments(
            preds = preds_df, 
            assignments = assign_df,
            trait = assign_traits[req_name], 
            assignment_col = assign_names[req_name],
            n_groups = 2, 
            jitter = (0.02, 0.01),  # jitter for (background, assigned) points; set None to turn off jitter
            gen = gen, 
            random_seed = 1,
            trait_name = None,
            outdir = req_outdir)
        
        density_assignments(
            preds = preds_df,
            assignments = assign_df,
            trait = assign_traits[req_name], 
            assignment_col = assign_names[req_name],
            n_groups = 2,
            gen = gen, 
            trait_name = None,
            outdir = req_outdir)    
    

def update_assignments(assignments, updates, all_requests, outdir=None):
    '''
    Update assignment files to reflect rats that were actually shipped to the assigned project.
    
    Args:
    
        assignments: A csv path or Pandas dataframe with previously proposed 
            assignment mapped to requests, as output by output_assignment_results().
        
        updates: Dict or path to an 'assignment_updates' json file. Must have 
            key 'requests' and values specifying paths to shipping sheets for 
            processed requests from the current generation.

        all_requests: Boolean. Whether or not all current requests are included
            in the update. If True, any RFIDs not found in the set of shipped 
            rats will be reclassified as 'not_assigned'. If False, updates will 
            only be made to the RFIDs included in the shipping sheets and 
            remaining assignments will remain unchanged.

        outdir: The directory in which to save output files.
                        
    Returns:
        A dataframe containing all updated assignment results.
    '''
    
    if isinstance(assignments, str):
        assign_df = pd.read_csv(assignments, dtype={'rfid': str})
    else:
        assign_df = assignments
    
    if isinstance(updates, str):
        with open(updates, 'r') as updates_file:
            updates = json.load(updates_file)
    else:
        updates = updates

    if not isinstance(updates, dict):
        raise TypeError("Input 'updates' must be a dictionary.")
    
    # make a datestamped output directory for all results
    if outdir is not None:
        if Path(outdir).parts[-1] != str(datestamp):
            outdir = os.path.join(outdir, datestamp)
        os.makedirs(outdir, exist_ok = True)

    # format dates
    assign_df['dob'] = pd.to_datetime(assign_df['dob'], format='mixed').dt.date
    assign_df['dow'] = pd.to_datetime(assign_df['dow'], format='mixed').dt.date

    # save the generation number
    gen = assign_df['generation'].tolist()[0]
    
    # organize output columns 
    assign_cols = ['rfid','animalid','earpunch','sex','coatcolor','generation',
                   'dob','dow','project_name','request_name','assignment', 'comments']
    col_dtypes = {'generation': 'Int64', 'animalid':'str', 'rfid':'str', 
                  'earpunch':'str', 'sex':'str', 'coatcolor':'str', 'dob':'str', 'dow':'str', 
                  'project_name':'str', 'request_name':'str', 'assignment':'str','comments':'str'}

    # names of requests to be updated
    requests = updates['requests']
    update_reqs = list(requests.keys())
    n_requests = len(update_reqs)
    assignments = []

    # list/set to check for double-assigned rfids
    double_assigned_rfids = []
    processed_rfids = set()

    for req_name in update_reqs:

        print(f'Processing request: {req_name} \n')
        # read in proposed assignments and shipped rats for the request
        req_dict = requests[req_name]
        req_type = req_dict['assignment_type']
        req_trt = req_dict['treatment']
        req_df = assign_df[assign_df['request_name'] == req_name]
        assigned_rfids = set(req_df['rfid'])
        proj_name = req_df['project_name'].tolist()[0]
        assignment = req_df['assignment'].tolist()[0]
        assignments.append(assignment)

        if req_type == 'hsw_breeders':
            pairs = pd.read_csv(req_dict['pairs_file'], \
                dtype={'dam_rfid':str,'sire_rfid':str})
            shipped_rfids = pairs['dam_rfid'].tolist() + \
                pairs['sire_rfid'].tolist()
        else:
            req_ss = os.path.join(updates['ss_dir'], \
                f'gen{updates['generation']}', req_dict['shipping_sheet'])
            req_ss = pd.read_excel(req_ss, dtype={'rfid': str})
            # subset to the assigned treatment if included in the shipping sheet
            if req_trt is not None:
                req_ss = req_ss[req_ss['treatment'] == req_trt]
            
            shipped_rfids = req_ss['rfid'].tolist()

        # add metadata into the request's dictionary
        req_dict['proj_name'] = proj_name
        req_dict['assignment'] = assignment
        req_dict['assign_df'] = req_df
        req_dict['assigned_rfids'] = assigned_rfids
        req_dict['shipped_rfids'] = shipped_rfids
        requests[req_name] = req_dict

        # check for duplicate assignments
        for rfid in shipped_rfids:
            if rfid in processed_rfids:
                double_assigned_rfids.append(rfid)
            else:
                processed_rfids.add(rfid)

    if len(double_assigned_rfids) > 0:
        print('Some RFIDs are double-assigned!')
        print(double_assigned_rfids)
        double_assigned_ids =  assign_df[assign_df['rfid']\
            .isin(double_assigned_rfids)]['animalid'].tolist()

        animalid_doubles = {animalid: None for animalid in double_assigned_ids}
        for animalid in double_assigned_ids:
            assigned_to_reqs = []
            for req in update_reqs:
                req_dict = requests[req]

                # get the shipped_rfids for THIS request (not the outer loop's)
                if req_dict['assignment_type'] == 'hsw_breeders':
                    pairs = pd.read_csv(req_dict['pairs_file'], 
                        dtype={'dam_rfid':str,'sire_rfid':str})
                    req_shipped_ids = pairs['dam_animalid'].tolist() + \
                        pairs['sire_animalid'].tolist()
                else:
                    req_ss_path = os.path.join(updates['ss_dir'], 
                        f'gen{updates["generation"]}', req_dict['shipping_sheet'])
                    req_ss_data = pd.read_excel(req_ss_path, dtype={'rfid': str})
                    req_trt_inner = req_dict['treatment']
                    if req_trt_inner is not None:
                        req_ss_data = req_ss_data[req_ss_data['treatment'] == req_trt_inner]
                    req_shipped_ids = req_ss_data['animalid'].tolist()

                if animalid in req_shipped_ids:
                    assigned_to_reqs.append(req)
            animalid_doubles[animalid] = assigned_to_reqs
        pprint(animalid_doubles)
        sys.exit()

    # reset all assignments to 'not_assigned'
    updated_df = assign_df.copy()
    updated_df['project_name'] = 'not_assigned'
    updated_df['request_name'] = 'not_assigned'
    updated_df['assignment'] = 'not_assigned'

    # mark assignments for each request
    for req_name in update_reqs:
        
        req_dict = requests[req_name]
        assignment = req_dict['assignment']
        proj_name = req_dict['proj_name']
        req_ids = req_dict['shipped_rfids'] # shipped IDs are the ACTUAL assignments

        req_idx = updated_df[updated_df['rfid'].isin(req_ids)].index

        updated_df.loc[req_idx, 'project_name'] = proj_name
        updated_df.loc[req_idx, 'request_name'] = req_name
        updated_df.loc[req_idx, 'assignment'] = assignment
        updated_df.loc[req_idx, assignment] = 1

        # save request-specific assignments to file
        req_df = updated_df.loc[req_idx,]
        req_df.sort_values('animalid', inplace=True)

        if outdir is not None:
            req_file = f'rattaca_gen{gen}_assign_{assignment}_{datestamp}.csv'
            req_file = os.path.join(outdir, req_file)
            req_df.to_csv(req_file, index=False, na_rep='')
            print(f'\n{assignment} assignments saved to {req_file}')

    # write all assignments to file
    updated_df.sort_values('animalid', inplace=True)

    assigned_rfids = updated_df[updated_df['assignment'] \
        =='not_assigned'][['rfid','animalid']]
    unassigned_rfids = updated_df[updated_df['assignment'] \
        != 'not_assigned'][['rfid','animalid']]

    if outdir is not None:
        updated_file = os.path.join(outdir, \
            f'rattaca_gen{gen}_assignments_{datestamp}.csv')
        assigned_rfids_file = \
            os.path.join(outdir, f'rattaca_gen{gen}_assigned_rfids_{datestamp}.csv')
        unassigned_rfids_file = \
            os.path.join(outdir, f'rattaca_gen{gen}_avail_rfids_{datestamp}.csv')
        
        updated_df.to_csv(updated_file, index=False, na_rep='')    
        print(f'\nAll updated assignments saved to {updated_file} \n')

        assigned_rfids.to_csv(assigned_rfids_file, index=False, na_rep='')
        unassigned_rfids.to_csv(unassigned_rfids_file, index=False, na_rep='')
    
    print('Final assignment counts:')
    print(updated_df[['assignment','sex']].value_counts())

    return(updated_df)


# UNTESTED
def permute_random():
    '''
    Assign rats to random projects using permutation to eliminate ordering biases.
    
    Args:
        all_requests: a list of random request objects to process.
        
    Returns:
        A dictionary mapping projects to assigned RFIDs.
    '''

    assignment_results = {}
    
    # continue until all requests are satisfied
    while any(not request.is_satisfied() for request in all_requests):
        # get unsatisfied requests
        open_requests = \
            [request for request in all_requests if not request.is_satisfied()]
        n_open_requests = len(open_requests)
                
        # try all possible permutations of request orders
        for permuted_order in permutations(range(n_open_requests)):

            # for each project request in the current permutation... 
            for current_request_idx in permuted_order:
                current_request = open_requests[current_request_idx]
                other_requests = [proj for i, proj in enumerate(open_requests) \
                    if i != current_request_idx]                
                
                # ...randomly choose one available rat for assignment
                assign_rfid = random.choice(current_request.available_rfids)

            if not current_request.is_satisfied():
                current_request.assign(assign_rfid, remaining_requests = other_requests)
                for rm_from_request in other_requests:
                    rm_from_request.remove(assign_rfid)
                
                # store results for this project
                project_name = current_request.project
                if project_name not in assignment_results:
                    assignment_results[project_name] = []
                assignment_results[project_name].extend(assign_rfids)
    
    print('all random assignments')
    print(assignment_results)
    return all_requests
