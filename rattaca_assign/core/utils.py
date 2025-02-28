"""
Utility functions for RATTACA assignment.
"""

import json
import pandas as pd
from rattaca_assign.models.request_rattaca import RATTACA
from rattaca_assign.models.request_random import RandProj
from rattaca_assign.models.request_breeders import HSWBreeders


# function to list which types of requests have been made
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
    for proj_request in args.requests:
        request = Request(proj_request, args)
        req_type = request.assignment_type
        which_types.append(req_type)
    return(list(set(which_types)))  



# function to clean up and format the colony dataframe
def prep_colony_df(args):
    '''
    Drops unusable samples from the colony dataframe.
    
    Drops all dead rats, rats with undetermined sex, rats lacking RFIDs,
    and excluded RFIDs from consideration for assignment. Adds one column
    denoting which rats have been genotyped (for use in assigning to 
    random-choice projects).

    Args:
        args.colony_dataframe: The path to the colony dataframe (csv).
        args.predictions: The path to the predictions csv.
        args.exclude: The path to the exclude list (csv).

    Returns:
        A pandas dataframe with colony data for all rats available 
        for assignment, with an added columns denoting which samples 
        have been genotyped.
    '''

    df = pd.read_csv(args.colony_dataframe[0], 
                                dtype = {'rfid': str, 'accessid': int})
    
    if args.predictions:
        preds_df = pd.read_csv(args.predictions[0], dtype = {'rfid': str})
        gtyped_rfids = preds_df['rfid'].tolist()
    else:
        gtyped_rfids = []

    # check for duplicate RFIDs
    if sum(df.duplicated(subset = ['rfid'])) > 0:
        print('Error: Colony metadata has duplicate RFIDs')
        print('Check the following samples before proceeding:')
        dups = df.duplicated(subset = ['rfid'])
        dups = df[dups]
        for rfid in dups['rfid']:
            print(rfid)
        exit()

    # drop rats without RFIDs
    df[~df['rfid'].isnull()]
    
    # identify dead rats
    dead_strs = ['dead', 'die', 'death', 'euth', 'eauth', 'kill']
    dead_search = '|'.join(dead_strs)
    dead_rats = df[df['comments']\
        .str.contains(dead_search, case=False, na=False)]['rfid'].tolist()
    
    # identify rats with unknown sex
    ambig_sex = df[df['comments']\
        .str.contains('unsure sex', case=False, na=False)]['rfid'].tolist()
    
    # drop dead rats, rats w/ unknown sex
    drop_rats = dead_rats + ambig_sex
    df = df[~df['rfid'].isin(drop_rats)]

    # identify rats that have been genotyped
    df['gtyped'] = df['rfid'].isin(gtyped_rfids).astype(int)

    return(df)


# function to assign breeders a priori
def prioritize_breeders(args):
    '''
    Assigns singleton offspring to HSW breeders by priority 
    prior to project assignments.
    
    Identifies all rats that are the only M or F from their litter,
    assigns them automatically to HSW breeders to ensure all breeder
    pairs have one M and one F represented in the next HSW generation.

    Args:
        args.colony_dataframe: The path to the colony dataframe (csv).
    
    Returns:
        A list of RFIDs that must be assigned to HSW breeders.
    '''
    
    df = prep_colony_df(args)
    m_breeders = []
    f_breeders = []

    for pair in df['breederpair'].unique():

        litter_males = \
            df[(df['breederpair'] == pair) & (df['sex'] == 'M')]['rfid']
        litter_females = \
            df[(df['breederpair'] == pair) & (df['sex'] == 'F')]['rfid']
        
        # any rat that is the single M or F from its litter 
        # is automatically assigned to HS West breeders
        if len(litter_males) == 1:
            assign_male = litter_males
            m_breeders.append(assign_male)
        
        if len(litter_females) == 1:
            assign_female = litter_females
            f_breeders.append(assign_female)
        
    breeders = m_breeders + f_breeders
    return(breeders)

def load_request_files(args):
    """
    Load request objects from JSON files.
    
    Args:
        args.requests: Command line arguments containing request files.
        
    Returns:
        A dictionary with separate lists of request objects for each 
        request type.
    """
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

