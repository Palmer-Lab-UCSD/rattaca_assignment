#!/usr/bin/env python3

### Use this script to assign HS West rats to RATTACA projects 
### Required inputs:
### 1. a list of RFIDs available for assignment
### 2. breeding values mapped to RFIDs for each trait/project being requested
### 3. sex values for each RFID
### 4. project requests specifying desired sample sizes, in .json format

### usage: python3 rattaca_assignment.py -c colony_dataframe.csv 
#       -p predictions.csv -o /output/path/ -r request1.json request2.json...


### dependencies ###
import json
import sys
import os
import pandas as pd
import numpy as np
import argparse
from itertools import permutations

### functions ###

def parse_args(args):
    parser = argparse.ArgumentParser(
                    prog='RATTACA assignment',
                    description='''Uses permutation to maximize per-project
                        between-group differences in predicted traits when
                        assigning rats to new projects''',
                    epilog='--PUT EPILOG HERE--')
    parser.add_argument('-o', '--output_prefix', 
        required=False, type=str, nargs=1, 
        dest="outfile_prefix", help="<output file prefix>")
    parser.add_argument('-c', '--colony_df', 
        required=True, type=str, nargs=1, 
        dest="colony_dataframe", help="<HS West colony dataframe csv>")
    parser.add_argument('-p', '--predictions', 
        required=False, type=str, nargs=1, 
        dest="predictions", help="<rattaca predictions csv")
    parser.add_argument('-r', '--request_files', 
        required=True, type=str, 
        nargs='+', dest="requests", help="<requests jsons")
    parser.add_argument('-e', '--exclude_rfids', 
        required=False, type=str, nargs=1,
        dest='exclude', help='RFIDs in colony data to exclude from assignment')

    return parser.parse_args(args)

# function to clean up and format the colony dataframe
def prep_colony_df(args):
    '''Drops unusable samples from the colony dataframe.
    
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
    '''Assigns singleton offspring to HSW breeders by priority 
    prior to project assignments.
    
    Identifies all rats that are the only M or F from their litter,
    assigns them automatically to HSW breeders to ensure all breeder
    pairs have one M and one F represented in the next HSW generation.

    Args:
        args.colony_dataframe: The path to the colony dataframe (csv).
        args.predictions: The path to the predictions csv.
    
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

# function to format results from all assignments to the database
def output_df_for_db(request_list):
    '''Extracts assignment data from multiple Request variables into a 
    dataframe formatted to upload to the Palmer Lab database.
    '''
    
    pass


# primary assignment algorithm
def main(args):
    
    # empty lists to hold different request types
    all_requests = []
    rattaca_requests = []
    random_requests = []

    # read in request files by assignment type
    for proj_request in args.requests:

        request = Request(proj_request, args)
        if request.assignment_type == 'rattaca':
            request = RATTACA(proj_request, args)
            rattaca_requests.append(request)
        if request.assignment_type == 'random':
            request = RandProj(proj_request, args)
            random_requests.append(request)
        if request.assignment_type == 'hsw_breeders':
            request = HSWBreeders(proj_request, args)
            hsw_breeders = request

    # assign HSW breeders by priority
    priority_breeders = prioritize_breeders(args)
    hsw_breeders.assign(priority_breeders)
    for request in [rattaca_requests + random_requests]:
        request.remove(priority_breeders)
        
    # list all requests for assignment
    all_requests = hsw_breeder_requests + rattaca_requests + random_requests

    # permute the assignment order to all requests
    ar = 0
    while any(not request.is_satisfied() for request in all_requests):

        open_requests = \
            [request for request in all_requests if not request.is_satisfied()]
        n_open_requests = len(open_requests)
        
        # initialize variables to store stats from the eventual best permutation
        best_delta = 0
        best_permutation = None
        best_rfids = None

        ar += 1
        print(f'ASSIGNMENT ROUND {ar}')
        for permuted_order in permutations(range(n_open_requests)):
            # print(f'permuted_order: {permuted_order}')
            proposed_rfids = {}
            permutation_delta = 0
            project_deltas = []

            # for each project request in the current permutation... 
            for current_request_idx in permuted_order:

                current_request = open_requests[current_request_idx]
                # print(f'current_request: {current_request.project}')
                # ...propose animals for assignment to the project 
                # and add the project's delta to the total delta for the permutation
                project_delta, proposal_rfids = current_request\
                    .proposal(proposed_rfids) 
                permutation_delta += project_delta
                project_deltas.append(project_delta)
                proposed_rfids[current_request_idx] = proposal_rfids

            # keep the proposed rfids for assignment from the permutation 
            # that maximizes delta
            if permutation_delta > best_delta:
                best_delta = permutation_delta
                best_permutation = permuted_order
                best_rfids = proposed_rfids

        # print(f'final best_permutation: {best_permutation}')
        # print(f'best_rfids: {best_rfids}')

        # assign rats to each request in the permuted order that maximizes delta
        for project_index in best_permutation:

            project_to_assign = open_requests[project_index]
            other_projects = [proj for i, proj in enumerate(open_requests) \
                              if i != project_index]
            
            assign_rfids = best_rfids[project_index]
            remove_rfids = [rfid for remaining_rfids in \
                            [rfid_val for idx_key, rfid_val in \
                             best_rfids.items() if idx_key != project_index] \
                                for rfid in remaining_rfids]
            print(f'assign_rfids: {assign_rfids}')
            print(f'remove_rfids: {remove_rfids}')

            if not project_to_assign.is_satisfied():
                project_to_assign.assign(assign_rfids)
                for rm_from_project in other_projects:
                    rm_from_project.remove(assign_rfids)
        # if ar == 3:
        #     break


### classes ###

class Request:
    
    def __init__(self, request_file, args):

        print('Initializing Request') # debug message - can delete later

        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadat = json.load(rf)
            self.assignment_type = req_metadat['assignment_type']
            self.project = req_metadat['project']
            self.trait = req_metadat['trait']
            print(f'Project: {self.project}, trait: {self.trait}') # debug message - delete later
            self.n_requested_total = req_metadat['n_rats']['total']
            self.n_requested_males = req_metadat['n_rats']['male']['total']
            self.n_requested_females = req_metadat['n_rats']['female']['total']
                
        # read in colony data 
        self.colony_df = prep_colony_df(args)

        # add 'trait' info to dataframe
        self.trait_metadata = prep_colony_df(args)
        self.trait_metadata[[self.trait]] = 0

        # create a list of rats available for assignment, ordered by trait prediction
        self.available_rfids = self.trait_metadata['rfid'].tolist()

        print('Finished initializing Request') # debug message - can delete later

    # function to remove assigned rats from the list of available rats
    def remove(self, rats_to_remove):
        
        if isinstance(rats_to_remove, list):
            rfids_to_remove = rats_to_remove
        elif isinstance(rats_to_remove, dict):
            rfids_to_remove = list(rats_to_remove.keys())
        else:
            raise TypeError('rats_to_remove must be either a list or a dictionary')

        for rfid in rfids_to_remove:
            
            if rfid in self.available_rfids:
                self.available_rfids.remove(rfid)

    # generic function to assign rats to a project, remove them from availability
    def assign(self, rats_to_assign):

        if self.assignment_type == 'rattaca':
            return self.assign_rattaca(rats_to_assign)
        elif self.assignment_type == 'random':
            return self.assign_random(rats_to_assign)
        elif self.assignment_type == 'hsw_breeders':
            return self.assign_hsw_breeders(rats_to_assign)
        else:
            raise NameError('Assignment type must be either "rattaca", \
                            "hsw_breeders", or "random". Please check the .json \
                            file for this project request')


    # generic function to check if the request has been fulfilled
    def is_satisfied(self, sex=None, group=None):
        
        if self.assignment_type == 'rattaca':
            return self.is_satisfied_rattaca(sex, group)
        elif self.assignment_type == 'random':
            return self.is_satisfied_random(sex, group)
        elif self.assignment_type == 'hsw_breeders':
            return self.is_satisfied_hsw_breeders(sex, group)
        else:
            raise NameError('''Assignment type must be either "rattaca", 
                            "hsw_breeders", or "random". Please check the .json 
                            file for this project request''')
            

### TO DO: flesh out sub-classes ##
# Request subclass for HS West breeders
class HSWBreeders(Request):
    
    def assign(self, rfids_to_assign):
        pass

    def is_satisfied_hsw_breeers(self):
        pass

# Request subclass for RATTACA projects
class RATTACA(Request):
    
    def __init__(self, request_file, args):
        
        print('Initializing RATTACA request') # debug message - delete later

        # inherit attributes from the parent Request class
        super().__init__(request_file, args)
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadat = json.load(rf)
            self.n_requested_males_high = req_metadat['n_rats']['male']['high']
            self.n_requested_males_low = req_metadat['n_rats']['male']['low']
            self.n_requested_females_high = \
                req_metadat['n_rats']['female']['high']
            self.n_requested_females_low = \
                req_metadat['n_rats']['female']['low']
            self.n_requested_high = \
                self.n_requested_males_high + self.n_requested_females_high
            self.n_requested_low = \
                self.n_requested_males_low + self.n_requested_females_low

        # initialize dictionaries to hold assigned rats
        self.assigned_males_high = {}
        self.assigned_males_low = {}
        self.assigned_females_high = {}
        self.assigned_females_low = {}
                
        # initialize the delta for the request
        self.delta = 0

        # merge trait predictions into colony data
        self.setup_trait_metadata(args)

        # create a list of rats available for assignment, ordered by trait prediction
        # self.available_rfids = self.trait_metadata['rfid'].tolist()

        print('Finished initializing RATTACA request') # debug message - delete later

    # function to set up trait metadata for a project request
    def setup_trait_metadata(self, args):
 
        if args.exclude:
            exclude_rfids = pd.read_csv(args.exclude[0], header = None)\
                .astype(str)[0].values
        else:
            exclude_rfids = []

        # read in predictions
        predictions_df = pd.read_csv(args.predictions[0], dtype = {'rfid': str})

        # merge predictions with colony metadata
        all_metadata = self.colony_df.merge(
            predictions_df[['rfid', self.trait]], on='rfid', how='right')

        # sort by the requested trait, remove NAs
        self.trait_metadata = all_metadata.sort_values(by = self.trait, 
            axis=0, ascending=False)\
                .dropna(subset=self.trait, ignore_index=True)

        # add a column of high/low group assignments
        self.trait_metadata[f'{self.trait}_2group'] = \
            self.trait_groups(self.trait, n_groups=2)

        # drop rats from the exclude list
        self.trait_metadata = self.trait_metadata[~self.trait_metadata['rfid'].\
            isin(exclude_rfids)]

        # drop any groups that aren't requested
        if self.n_requested_males == 0:
            self.trait_metadata = self.trait_metadata\
                [self.trait_metadata['sex'] != 'M']
        if self.n_requested_females == 0:
            self.trait_metadata = self.trait_metadata\
                [self.trait_metadata['sex'] != 'F']
        if self.n_requested_high == 0:
            self.trait_metadata = self.trait_metadata\
                [self.trait_metadata[f'{self.trait}_2group'] != 'high']
        if self.n_requested_low == 0:
            self.trait_metadata = self.trait_metadata\
                [self.trait_metadata[f'{self.trait}_2group'] != 'low']

        # create a list of rats available for assignment, ordered by trait prediction
        self.available_rfids = self.trait_metadata['rfid'].tolist()


    # function to propose an assignment of RFIDs to a project
    # unavail_rats = either empty, or a list of RFIDs,
    # or a dict of (sex, pred) w/ RFIDs as keys
    def proposal(self, unavail_rats = None):
        '''Identifies the two most extreme available rats and calculates
        the difference (delta) between their predictions.
        
        Samples the two rats from the list of available RFIDs with the
        highest and lowest predicted trait value, which become the rats
        'proposed' for assignment to the request. Calculates the delta 
        between these predictions and adds this to the actual total delta
        for the request to calculate the proposed delta for the request.

        Args:
            unavail_rats: An optional list of RFIDs or dictionary with 
            RFID keys denoting which rats should be excluded from the 
            proposed assignment.
            
            self.available_rfids: The list of RFIDs currently available
            for assignment to the request.

        Returns:
            A tuple with two elements: 
                [0]: The proposed delta for the proposed assignment
                [1]: The list of RFIDs of the two proposed rats to assign, 
                sorted in order [lowest prediction, highest prediction]
        '''
        if unavail_rats is None:
            unavail_rfids = []
        elif isinstance(unavail_rats, list):
            unavail_rfids = unavail_rats
        elif isinstance(unavail_rats, dict):
            unavail_rfids = list(unavail_rats.keys())
        else:
            raise TypeError('unavail_rats must be either a list, a dictionary, or None')
        
        # initial indexes to start search for extreme available rats
        i = 0    # first element is the rat with the max trait value
        j = len(self.available_rfids) - 1 # last element has the min trait value

        # select the most extreme high & low available rats
        while self.available_rfids[i] in unavail_rfids:
            i += 1
        while self.available_rfids[j] in unavail_rfids:
            j -= 1
        max_rat = self.available_rfids[i]
        min_rat = self.available_rfids[j]

        # calculate delta for this iteration: the difference 
        # in predicted trait values
        min_rat_sex, min_rat_pred, min_rat_group = self.available_rats[min_rat]
        max_rat_sex, max_rat_pred, max_rat_group = self.available_rats[max_rat]
        delta_iter = max_rat_pred - min_rat_pred

        # calculate the proposed total delta: the sum of all iterations' 
        # deltas if this proposed iteration were to be added
        proposed_delta = self.delta + delta_iter

        # return proposed delta, selected rats
        return (proposed_delta, [min_rat, max_rat])

    # function to assign rats to RATTACA projects based on their genetic predictions
    def assign_rattaca(self, rfids_to_assign): # rats_to_assign = best_rfids

        # get predictions directly from metadata 
        # (instead of available_rats - will give an uninformative key error)
        rfid_low = rfids_to_assign[0]
        rfid_high = rfids_to_assign[1]
        pred_low = self.trait_metadata\
            [self.trait_metadata['rfid'] == rfid_low][self.trait].iloc[0]
        pred_high = self.trait_metadata\
            [self.trait_metadata['rfid'] == rfid_high][self.trait].iloc[0]

        if (isinstance(rfids_to_assign, list)) & (pred_low < pred_high):
            pass
        else:
            raise TypeError('''rfids_to_assign must be a list with two RFID ''' + 
            '''elements in order [min_rat, max_rat], as output by proposal()''')
        
        rats_to_assign = {}
        # ensure rats designated for assignment are actually available,
        # get relevant metadata for available RFIDS
        ### TO DO: test whether it's better to return None or raise value error here
        for rfid in rfids_to_assign:
            if rfid in self.available_rfids:
                rats_to_assign[rfid] = self.available_rats[rfid]
            else:
                print(f'''RFID {rfid} is not available for ''' +
                '''assignment, was NOT assigned''')
                return None

        min_rat = rfids_to_assign[0]
        max_rat = rfids_to_assign[1]
        min_rat_sex, min_rat_pred, min_rat_group = rats_to_assign[min_rat]
        max_rat_sex, max_rat_pred, max_rat_group = rats_to_assign[max_rat]

        # assign the high rat and remove it from availability
        if max_rat_sex == 'M':
            assign_to = self.assigned_males_high
        elif max_rat_sex == 'F':
            assign_to = self.assigned_females_high
        
        # if the request is not yet satisfied
        # assign the high rat and remove it from availability
        if self.is_satisfied(max_rat_sex, max_rat_group):
            message = (f'RFID {max_rat} could not be assigned: '
                       f'All {max_rat_sex}/{max_rat_group} group '
                       'assignments already satisfied')
            print(message)
        else:
            assign_to[max_rat] = rats_to_assign[max_rat]
            self.remove([max_rat])
                
        # set up assignment for the low rat and remove it from availability
        if min_rat_sex == 'M':
            assign_to = self.assigned_males_low
        elif min_rat_sex == 'F':
            assign_to = self.assigned_females_low

        # if the request is not yet satisfied
        # assign the low rat and remove it from availability
        if self.is_satisfied(min_rat_sex, min_rat_group):
            message = (f'RFID {min_rat} could not be assigned: '
            f'All {min_rat_sex}/{min_rat_group} group '
            'assignments already satisfied')
            print(message)
        else:
            assign_to[min_rat] = rats_to_assign[min_rat]
            self.remove([min_rat])
        
            # get the delta betwen assigned rats to update the request's total delta
            min_rat_pred = rats_to_assign[min_rat][1]
            max_rat_pred = rats_to_assign[max_rat][1]
            assignment_delta = max_rat_pred - min_rat_pred
            self.delta += assignment_delta

        # automatically update list of available rats once a given subsample
        # for a request has been satisfied
        if self.is_satisfied('M','high'):
            drop_df = self.trait_metadata
            drop_rfids = drop_df[(drop_df['sex']=='M') & \
                        (drop_df[f'{self.trait}_2group']=='high')]['rfid'].tolist()
            self.available_rfids = [rfid for rfid in self.available_rfids \
                                    if rfid not in drop_rfids]
        if self.is_satisfied('M','low'):
            drop_df = self.trait_metadata
            drop_rfids = drop_df[(drop_df['sex']=='M') & \
                        (drop_df[f'{self.trait}_2group']=='low')]['rfid'].tolist()
            self.available_rfids = [rfid for rfid in self.available_rfids \
                                    if rfid not in drop_rfids]
        if self.is_satisfied('F','high'):
            drop_df = self.trait_metadata
            drop_rfids = drop_df[(drop_df['sex']=='F') & \
                        (drop_df[f'{self.trait}_2group']=='high')]['rfid'].tolist()
            self.available_rfids = [rfid for rfid in self.available_rfids \
                                    if rfid not in drop_rfids]
        if self.is_satisfied('F','low'):
            drop_df = self.trait_metadata
            drop_rfids = drop_df[(drop_df['sex']=='F') & \
                        (drop_df[f'{self.trait}_2group']=='low')]['rfid'].tolist()
            self.available_rfids = [rfid for rfid in self.available_rfids \
                                    if rfid not in drop_rfids]


    # function to manually assign rats to a projects as needed
    def assign_manual(self, rfids_to_assign): # rats_to_assign = best_rfids

        if isinstance(rfids_to_assign, list):
            pass
        else:
            raise TypeError('rfids_to_assign must be a list')
        
        rats_to_assign = {}
        # ensure rats designated for assignment are actually available,
        # get relevant metadata for available RFIDS
        for rfid in rfids_to_assign:
            if rfid not in self.available_rfids:
                raise ValueError(f'''RFID {rfid} is not available for assignment, 
                                 was NOT assigned''')
            else:
                rats_to_assign[rfid] = self.available_rats[rfid]        
            rat_sex =  rats_to_assign[rfid][0]
            rat_group =  rats_to_assign[rfid][2]

            if (rat_sex == 'M') & (rat_group == 'high'):
                assign_to = self.assigned_males_high
                n_requested = self.n_requested_males_high
            elif (rat_sex == 'M') & (rat_group == 'low'):
                assign_to = self.assigned_males_low
                n_requested = self.n_requested_males_low
            elif (rat_sex == 'F') & (rat_group == 'high'):
                assign_to = self.assigned_females_high
                n_requested = self.n_requested_females_high
            elif (rat_sex == 'F') & (rat_group == 'low'):
                assign_to = self.assigned_females_low
                n_requested = self.n_requested_females_low

            if not self.is_satisfied(rat_sex, rat_group):
                assign_to[rfid] = rats_to_assign[rfid]
                self.remove([rfid])
            else:
                print(f'''RFID {rfid} could not be assigned. All \
                      {n_requested} {assign_to} assignments are \
                        already satisfied''')


    # function to check if a prediction-based RATTACA project has 
    # successfully assigned all requested rats
    def is_satisfied_rattaca(self, sex=None, group=None):

        if (sex is None) & (group is None):
            n_remaining = self.n_remaining['n_remaining_total']
        if (sex == 'M') & (group is None):
            n_remaining = self.n_remaining['n_remaining_males']
        elif (sex == 'F') & (group is None):
            n_remaining = self.n_remaining['n_remaining_females']
        if (group == 'high') & (sex is None):
            n_remaining = self.n_remaining['n_remaining_high']
        if (group == 'low') & (sex is None):
            n_remaining = self.n_remaining['n_remaining_low']
        if (sex == 'M') & (group == 'high'):
            n_remaining = self.n_remaining['n_remaining_males_high']
        elif (sex == 'M') & (group == 'low'):
            n_remaining = self.n_remaining['n_remaining_males_low']
        elif (sex == 'F') & (group == 'high'):
            n_remaining = self.n_remaining['n_remaining_females_high']
        elif (sex == 'F') & (group == 'low'):
             n_remaining = self.n_remaining['n_remaining_females_low']

        return(n_remaining == 0)

    # function to assign high/low group values
    def trait_groups(self, trait, n_groups=2):

        df = self.trait_metadata
        preds = df[trait].tolist()
        trait_quantiles = np.quantile(a = preds, 
                                      q = np.linspace(0,1, n_groups+1))
        if n_groups == 2:
            def get_2group(val, quantiles):
                if val >= quantiles[1]:
                    return 'high'
                else:
                    return 'low'
            groups = [get_2group(pred, trait_quantiles) for pred in preds]
        
        if n_groups == 3:
            def get_3group(val, quantiles):
                if val < quantiles[1]:
                    return 'low'
                elif val >= quantiles[2]:
                    return 'high'
                else:
                    return 'mid'
            groups = [get_3group(pred, trait_quantiles) for pred in preds]
        
        if n_groups > 3:
            groups = np.digitize(preds, trait_quantiles)
            
        return groups

    # property to keep track of all rats currently available for assignment
    @property
    def available_rats(self):

        # update trait metadata to include only rats that are currently 
        # available for assignment
        current_metadata = self.trait_metadata\
            [self.trait_metadata['rfid'].isin(self.available_rfids)]
        # convert the metadata df to a dictionary 
        # with key:RFID and value:(sex, prediction, group)
        current_metadata = current_metadata.\
            set_index('rfid')[['sex', self.trait, f'{self.trait}_2group']].\
            to_dict(orient='index')
        available_rats = \
            {k: (v['sex'], v[self.trait], v[f'{self.trait}_2group']) \
                for k, v in current_metadata.items()}
        
        return available_rats
    
    # property to track all males currently available for assignment
    @property
    def available_males(self):
        return {k: v for k, v in self.available_rats.items() if v[0] == 'M'}
    
    # property to track all females currently available for assignment
    @property
    def available_females(self):
        return {k: v for k, v in self.available_rats.items() if v[0] == 'F'}

    # properties to track rats that have been assigned to the project
    @property
    def assigned_rats(self):
        assigned_high = self.assigned_males_high | self.assigned_females_high
        assigned_low = self.assigned_males_low | self.assigned_females_low
        assigned_males = self.assigned_males_high | self.assigned_males_low
        assigned_females = self.assigned_females_high | self.assigned_females_low
        assigned_total = assigned_males | assigned_females
        out = {'assigned_high': assigned_high, 
               'assigned_low': assigned_low,
               'assigned_males': assigned_males, 
               'assigned_females': assigned_females,
               'assigned_total': assigned_total}
        return out

   # counters to track the number of assigned rats
    @property
    def n_assigned(self):
        n_assigned_males_high = len(self.assigned_males_high)
        n_assigned_males_low = len(self.assigned_males_low)
        n_assigned_females_high = len(self.assigned_females_high)
        n_assigned_females_low = len(self.assigned_females_low)
        n_assigned_high = len(self.assigned_rats['assigned_high'])
        n_assigned_low = len(self.assigned_rats['assigned_low'])
        n_assigned_males = len(self.assigned_rats['assigned_males'])
        n_assigned_females = len(self.assigned_rats['assigned_females'])
        n_assigned_total = len(self.assigned_rats['assigned_total'])

        out = {'n_assigned_males_high': n_assigned_males_high,
               'n_assigned_males_low': n_assigned_males_low,
               'n_assigned_females_high': n_assigned_females_high,
               'n_assigned_females_low': n_assigned_females_low,
               'n_assigned_high': n_assigned_high,
               'n_assigned_low': n_assigned_low,
               'n_assigned_males': n_assigned_males,
               'n_assigned_females': n_assigned_females,
               'n_assigned_total': n_assigned_total}
        
        return out

    # counters to track the number of assignments remaining
    @property
    def n_remaining(self):

        n_remaining_males_high = \
            self.n_requested_males_high - self.n_assigned['n_assigned_males_high']
        n_remaining_males_low = \
            self.n_requested_males_low - self.n_assigned['n_assigned_males_low']
        n_remaining_females_high = \
            self.n_requested_females_high - self.n_assigned['n_assigned_females_high']
        n_remaining_females_low = \
            self.n_requested_females_low - self.n_assigned['n_assigned_females_low']
        n_remaining_high = self.n_requested_high - self.n_assigned['n_assigned_high']
        n_remaining_low = self.n_requested_low - self.n_assigned['n_assigned_low']
        n_remaining_males = self.n_requested_males - self.n_assigned['n_assigned_males']
        n_remaining_females = self.n_requested_females - self.n_assigned['n_assigned_females']
        n_remaining_total = self.n_requested_total - self.n_assigned['n_assigned_total']

        out = {'n_remaining_males_high': n_remaining_males_high,
               'n_remaining_males_low': n_remaining_males_low,
               'n_remaining_females_high': n_remaining_females_high,
               'n_remaining_females_low': n_remaining_females_low,
               'n_remaining_high': n_remaining_high,
               'n_remaining_low': n_remaining_low,
               'n_remaining_males': n_remaining_males,
               'n_remaining_females': n_remaining_females,
               'n_remaining_total': n_remaining_total}
        
        return out

# Request subclass for randomly-assigned projects
class RandProj(Request):
    
    def __init__(self, request_file, args, sibs_per_litter=None):
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadat = json.load(rf)
            self.assignment_type = req_metadat['assignment_type']
            self.project = req_metadat['project']
            self.trait = req_metadat['trait']
            self.n_requested_total = req_metadat['n_rats']['total']
            self.n_requested_males = req_metadat['n_rats']['male']['total']
            self.n_requested_females = req_metadat['n_rats']['female']['total']

        # initialize dictionaries to hold assigned rats
        self.assigned_males = {}
        self.assigned_females = {}

    def is_satisfied(self, n_sibs_allowed=0):
        pass


### execute ###

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)