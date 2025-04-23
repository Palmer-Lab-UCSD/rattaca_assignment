"""
Utility functions for RATTACA assignment that rely on request model classes.
"""

import json
import pandas as pd
from rattaca_assign.core.base_utils import prep_colony_df


def load_request_files(args):
    """
    Load request objects from JSON files.
    
    Args:
        args.requests: Command line arguments containing request files.
        
    Returns:
        A dictionary with separate lists of request objects for each 
        request type.
    """
    # import here to avoid circular imports
    from rattaca_assign.models.request_rattaca import RATTACA
    from rattaca_assign.models.request_random import RandProj
    from rattaca_assign.models.request_breeders import HSWBreeders
    
    rattaca_requests = []
    random_requests = []
    breeder_requests = []
    
    for proj_request in args.requests:
        with open(proj_request, 'r') as rf:
            req_metadata = json.load(rf)
            req_type = req_metadata['assignment_type']
            
            if req_type == 'rattaca':
                request = RATTACA(proj_request, args)
                rattaca_requests.append(request)
            elif req_type == 'random':
                request = RandProj(proj_request, args)
                random_requests.append(request)
            elif req_type == 'hsw_breeders':
                request = HSWBreeders(proj_request, args)
                breeder_requests.append(request)
    
    return {
        'rattaca': rattaca_requests,
        'random': random_requests,
        'hsw_breeders': breeder_requests
    }

def which_requests(args):
    '''
    List which types of requests have been made.
    
    Args:
        args.requests: Command line arguments containing paths to 
        request files.
        
    Returns:
        A list of unique request types found in the request files.
    '''
    which_types = []
    for request_file in args.requests:
        with open(request_file, 'r') as rf:
                req_metadata = json.load(rf)
                req_type = req_metadata['assignment_type']
                which_types.append(req_type)
    return(set(which_types))

# probably drop - superceded by HSWBreeders.update_priority_breeders()
# def _prioritize_breeders(args, requests):
#     '''
#     Assigns singleton offspring to HSW breeders by priority AFTER 
#     project assignments have been initiated.
    
#     Identifies all unassigned rats that are the only M or F from their litter,
#     assigns them automatically to HSW breeders to ensure all breeder
#     pairs have one M and one F represented in the next HSW generation. Singleton 
#     rats are identified by counting family members assigned to ALL current 
#     requests. This function must be run automatically after each assignment made 
#     to any non-breeder request.

#     Args:
#         args.colony_dataframe: The path to the colony dataframe (csv).
#         requests: A list of all request objects (both open and fulfilled) 
#             from the current generation.
    
#     Returns:
#         A list of RFIDs that must be assigned to HSW breeders.
#     '''
    
#     df = prep_colony_df(args)

#     for pair in df['breederpair'].unique():

#         litter_males = \
#             df[(df['breederpair'] == pair) & (df['sex'] == 'M')]['rfid']
#         litter_females = \
#             df[(df['breederpair'] == pair) & (df['sex'] == 'F')]['rfid']
        
#     # list all rats that have been assigned
#     assigned_males = []
#     assigned_females = []
#     for request in requests:
#         request_males = list(request.assigned_rats['assigned_males'].keys())
#         request_females = list(request.assigned_rats['assigned_females'].keys())
#         assigned_males.extend(request_males)
#         assigned_females.extend(request_females)

#     male_breeders = list(self.assigned_males.keys())
#     female_breeders = list(self.assigned_females.keys())


def prioritize_breeders(breeder_request, non_breeder_requests):
    '''
    Identify singleton offspring to assign to HSW breeders by priority AFTER 
    project assignments have been initiated.
    
    Identifies all unassigned rats that are the only M or F from their 
    litter, and returns them for assignment to HSW breeders to ensure all 
    current breeder pairs have one M and one F represented in the next HSW 
    generation. Singleton rats are identified by counting family members 
    assigned to ALL current requests. This function must be run 
    automatically after each assignment made to any non-breeder request.

    Args:
        args.colony_dataframe: The path to the colony dataframe (csv).
        breeder_request: An HSWBreeders request object
        non_breeder_requests: A list of all non-request objects (both open and 
            fulfilled) from the current generation.
    
    Returns:
        A list of RFIDs that must be assigned to HSW breeders.

    '''
    if isinstance(non_breeder_requests, list):
        pass
    else:
        raise TypeError('non_breeder_requests must be a list of Request classes')

    # list all rats that have been assigned to other projects
    assigned_males = []
    assigned_females = []
    for request in non_breeder_requests:
        request_males = \
            list(request.assigned_rats['assigned_males'].keys())
        request_females = \
            list(request.assigned_rats['assigned_females'].keys())
        assigned_males.extend(request_males)
        assigned_females.extend(request_females)

    male_breeders = list(breeder_request.assigned_males.keys())
    female_breeders = list(breeder_request.assigned_females.keys())
    
    # df = pd.read_csv(args.colony_dataframe[0], 
    #                             dtype = {'rfid': str, 'accessid': int})
    df = breeder_request.colony_df
    males_to_assign = []
    females_to_assign = []

    # count the number of unassigned rats remaining from each family
    for pair in df['breederpair'].unique():
        
        litter_males = \
            df[(df['breederpair'] == pair) & (df['sex'] == 'M')]['rfid']
        litter_females = \
            df[(df['breederpair'] == pair) & (df['sex'] == 'F')]['rfid']
        n_litter_males = len(litter_males)
        n_litter_females = len(litter_females)

        # if no males from the litter have yet been assigned as breeders...
        if not bool(set(male_breeders) & set(litter_males)):

            # ... count the number of available males remaining
            unassigned_males = set(litter_males) - set(assigned_males)

            # if only one male remains available, assign it to breeders
            if len(unassigned_males) == 1:
                males_to_assign.append(list(unassigned_males)[0])
        
        # check female assignments, same as males above
        if not bool(set(female_breeders) & set(litter_females)):
            unassigned_females = set(litter_females) - set(assigned_females)

            if len(unassigned_females) == 1:
                females_to_assign.append(list(unassigned_females)[0])

    priority_breeders = males_to_assign + females_to_assign

    return priority_breeders
