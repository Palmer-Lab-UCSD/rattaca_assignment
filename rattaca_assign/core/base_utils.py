"""
Basic utility functions that don't depend on request model classes.
"""

import pandas as pd

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
    dead_strs = ['dead', 'die', 'death', 'euth', 'eauth', 'kill', 'starve']
    dead_search = '|'.join(dead_strs)
    dead_rats = df[df['comments']\
        .str.contains(dead_search, case=False, na=False)]['rfid'].tolist()
    
    # identify rats with unknown sex
    ambig_sex = df[df['comments']\
        .str.contains('unsure sex', case=False, na=False)]['rfid'].tolist()
    
    # identify flooded cages
    flood_strs = ['flood', 'drown']
    flood_search = '|'.join(flood_strs)
    flooded_fams = df[df['comments']\
        .str.contains(flood_search, case=False, na=False)]['breederpair'].tolist()
    flooded_rats = df[df['breederpair'].isin(flooded_fams)]['rfid'].tolist()
    
    # drop dead rats, rats w/ unknown sex, rats from flooded cages
    drop_rats = dead_rats + ambig_sex + flooded_rats
    df = df[~df['rfid'].isin(drop_rats)]

    # identify rats that have been genotyped
    df['gtyped'] = df['rfid'].isin(gtyped_rfids).astype(int)

    return(df)


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
            df[(df['breederpair'] == pair) & \
                (df['sex'] == 'M')]['rfid'].tolist()
        litter_females = \
            df[(df['breederpair'] == pair) & \
                (df['sex'] == 'F')]['rfid'].tolist()
        
        # any rat that is the single M or F from its litter 
        # is automatically assigned to HS West breeders
        if len(litter_males) == 1:
            assign_male = litter_males[0]
            m_breeders.append(assign_male)
        
        if len(litter_females) == 1:
            assign_female = litter_females[0]
            f_breeders.append(assign_female)
    
    breeders = m_breeders + f_breeders
    return(breeders)