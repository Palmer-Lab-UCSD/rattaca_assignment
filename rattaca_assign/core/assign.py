'''
This is RATTACA's core assignment algorithm.

This module contains the main assignment logic that distributes rats to 
different projects based on their genetic predictions and other criteria.
'''

import pandas as pd
from random import choice
from itertools import permutations
from rattaca_assign.core.base_utils import prioritize_breeders
from rattaca_assign.core.model_utils import load_request_files, which_requests


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
    
    # load request objects
    request_objects = load_request_files(args)
    rattaca_requests = request_objects['rattaca']
    random_requests = request_objects['random']
    breeder_requests = request_objects['hsw_breeders']
    
    # STAGE 1: priority assignment to HSW breeders
    breeders_requested = 'hsw_breeders' in requests_made
    if breeders_requested:
        priority_breeders = prioritize_breeders(args)
        print(f'priority_breeders: {priority_breeders}')

        # assign breeders by priority
        for request in breeder_requests:
            request.assign_manual_hsw_breeders(priority_breeders)
        
        # remove assigned breeders from availability for other projects
        for request in rattaca_requests + random_requests:
            request.remove(priority_breeders)
    
    # STAGE 2: RATTACA assignments
    if rattaca_requests:
        # continue assigning RATTACA rats until all RATTACA requests are satisfied
        rattaca_assignment(rattaca_requests, breeder_requests)
    
    # STAGE 3: Complete breeder assignments after RATTACA
    if breeders_requested:
        for request in breeder_requests:
            if not request.is_satisfied():
                # Find any new priority breeders after RATTACA assignments
                new_priority_breeders = breeder_req.prioritize_breeders(args)
                if new_priority_breeders:
                    breeder_req.assign_manual_hsw_breeders(new_priority_breeders)
                    # Remove from random projects
                    for request in random_requests:
                        request.remove(new_priority_breeders)
                
                # For families that have contributed to RATTACA, ensure a male and female 
                # from each family goes to breeding if available
                complete_breeding_assignments(breeder_req)
    
    # STAGE 4: Only now assign to random projects
    if random_requests:
        # At this point, all RATTACA and priority breeder assignments are complete
        # Random projects get what's left
        for random_req in random_requests:
            # Ensure random projects only use rats from families that have already
            # contributed to RATTACA or breeders when possible
            assign_to_random_projects(random_req, rattaca_requests, breeder_requests)
    
    return rattaca_requests + breeder_requests + random_requests


def rattaca_assignment(rattaca_requests, non_rattaca_requests = None):
    '''
    Run RATTACA assignment algorithm until all requests are satisfied.
    
    Args:
        rattaca_requests: List of RATTACA request objects
        breeder_requests: List of breeder request objects (for updating)
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


def output_assignment_results(assignment_results, args, output_prefix=None):
    '''
    Output assignment results to CSV files.
    
    Args:
        assignment_results: a dictionary mapping projects to assigned RFIDs.
        output_prefix: path and file basename prefix to write output files.
        
    Returns:
        A dataframe containing all assignment results.
    '''

    # read in the colony dataframe
    # df = pd.read_csv(args.colony_dataframe[0], 
    #                             dtype = {'rfid': str, 'accessid': int})
                                
    # # combine all results into a single dataframe
    # all_assignments = []
    # for project, rfids in assignment_results.items():
    #     for rfid in rfids:
    #         all_assignments.append({
    #             'project': project,
    #             'rfid': rfid
    #         })
    
    for req in assignment_results:
        print('req.assigned_rats:')
        print(list(req.assigned_rats['assigned_total'].keys()))

    # df = pd.DataFrame(all_assignments)
    
    # write all assignments to file
    if output_prefix:
        assignments_file = f'{output_prefix}_all_assignments.csv'
        df.to_csv(assignments_file, index=False)
        
        # create per-project files
        print(f'\n')
        for project in df['project'].unique():
            project_df = df[df['project'] == project]
            project_file = f'{output_prefix}_{project}_assignments.csv'
            project_df.to_csv(project_file, index=False)
            print(f'{project} assignments saved to {project_file}')
        
        print(f'All assignments saved to {output_prefix}_all_assignments.csv \n')
    
    return df

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
