'''
This is RATTACA's core assignment algorithm.

This module contains the main assignment logic that distributes rats to 
different projects based on their genetic predictions and other criteria.
'''

import pandas as pd
from random import choice
from itertools import permutations
from rattaca_assign.core.model_utils import * 

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
    df = pd.read_csv(args.colony_dataframe[0], 
                     dtype={'rfid': str, 'accessid': int})
                                
    all_assignments = []
    
    # process each request object and extract assigned rats
    for req in assigned_requests:
        project_name = req.project
        
        # add assigned rats to the results
        assigned_rfids = list(req.assigned_rats.keys())
        for rfid in assigned_rfids:
            all_assignments.append({
                'project': project_name,
                'rfid': rfid
            })
    
    # create dataframe from all assignments
    result_df = pd.DataFrame(all_assignments)
    
    # write all assignments to file
    if output_prefix:
        assignments_file = f'{output_prefix}_all_assignments.csv'
        result_df.to_csv(assignments_file, index=False)
        
        # create per-project files
        print(f'\n')
        for project in result_df['project'].unique():
            project_df = result_df[result_df['project'] == project]
            project_file = f'{output_prefix}_{project}_assignments.csv'
            project_df.to_csv(project_file, index=False)
            print(f'{project} assignments saved to {project_file}')
        
        print(f'All assignments saved to {assignments_file} \n')
    
    return result_df


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
