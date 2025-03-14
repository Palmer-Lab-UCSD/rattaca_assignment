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

# # Move prioritize_breeders without model dependencies
# def prioritize_breeders(args):
#     '''
#     Assigns singleton offspring to HSW breeders by priority 
#     prior to project assignments.
    
#     Identifies all rats that are the only M or F from their litter,
#     assigns them automatically to HSW breeders to ensure all breeder
#     pairs have one M and one F represented in the next HSW generation.

#     Args:
#         args.colony_dataframe: The path to the colony dataframe (csv).
    
#     Returns:
#         A list of RFIDs that must be assigned to HSW breeders.
#     '''
    
#     df = prep_colony_df(args)
#     m_breeders = []
#     f_breeders = []

#     for pair in df['breederpair'].unique():

#         litter_males = \
#             df[(df['breederpair'] == pair) & (df['sex'] == 'M')]['rfid']
#         litter_females = \
#             df[(df['breederpair'] == pair) & (df['sex'] == 'F')]['rfid']
        
#         # any rat that is the single M or F from its litter 
#         # is automatically assigned to HS West breeders
#         if len(litter_males) == 1:
#             assign_male = litter_males
#             m_breeders.append(assign_male)
        
#         if len(litter_females) == 1:
#             assign_female = litter_females
#             f_breeders.append(assign_female)
        
#     breeders = m_breeders + f_breeders
#     return(breeders)