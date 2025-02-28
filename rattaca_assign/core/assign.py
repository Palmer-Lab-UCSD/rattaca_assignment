'''
This is RATTACA's core assignment algorithm.

This module contains the main assignment logic that distributes rats to 
different projects based on their genetic predictions and other criteria.
'''

import pandas as pd
from itertools import permutations
from rattaca_assign.core.utils import prioritize_breeders, load_request_objects, which_requests


def run_assignments(args):
    '''
    Run the assignment algorithm to distribute rats to projects.
    
    Args:
        args: command line arguments containing request files.
        
    Returns:
        A dictionary mapping projects to assigned RFIDs.
    '''

    # get request types
    requests_made = which_requests(args)
    
    # load request objects
    request_objects = load_request_files(args)
    rattaca_requests = request_objects['rattaca']
    random_requests = request_objects['random']
    breeder_requests = request_objects['hsw_breeders']
    
    # build list of all requests to process
    all_requests = []
    if 'rattaca' in requests_made:
        all_requests.extend(rattaca_requests)
    if 'random' in requests_made:
        all_requests.extend(random_requests)
    
    # handle HSW breeders if requested
    if 'hsw_breeders' in requests_made:

        # assign HSW breeders by priority
        priority_breeders = prioritize_breeders(args)
        for breeder_req in breeder_requests:
            breeder_req.assign(priority_breeders)
            
        # remove priority breeders from consideration for other projects
        for request in all_requests:
            request.remove(priority_breeders)
            
        # add breeder requests to all requests
        all_requests.extend(breeder_requests)
    
    # run the assignment algorithm
    return permute_assignment(all_requests)


def permute_assignment(all_requests):
    '''
    Assign rats to RATTACA projects using permutation to maximize overall delta.
    
    Args:
        all_requests: a list of request objects to process.
        
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
    
    return assignment_results


def output_assignment_results(assignment_results, output_path=None):
    '''
    Output assignment results to CSV files.
    
    Args:
        assignment_results: a dictionary mapping projects to assigned RFIDs.
        output_path: path to write output files.
        
    Returns:
        A dataframe containing all assignment results.
    '''
    # combine all results into a single dataframe
    all_assignments = []
    for project, rfids in assignment_results.items():
        for rfid in rfids:
            all_assignments.append({
                'project': project,
                'rfid': rfid
            })
    
    df = pd.DataFrame(all_assignments)
    
    # write all assignments to file
    if output_path:
        df.to_csv(f'{output_path}/all_assignments.csv', index=False)
        
        # create per-project files
        for project in df['project'].unique():
            project_df = df[df['project'] == project]
            project_df.to_csv(f'{output_path}/{project}_assignments.csv', \
                index=False)
    
    return df