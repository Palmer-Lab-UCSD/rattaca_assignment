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
        required=True, type=str, nargs=1, 
        dest="predictions", help="<rattaca predictions csv")
    parser.add_argument('-r', '--request_files', 
        required=True, type=str, 
        nargs='+', dest="requests", help="<requests jsons")
    parser.add_argument('-e', '--exclude_rfids', 
        required=False, type=str, nargs=1,
        dest='exclude', help='RFIDs in colony data to exclude from assignment')

    return parser.parse_args(args)


def main(args):
    
    # empty lists to hold different request types
    all_requests = []
    rattaca_requests = []
    hsw_breeder_requests = []
    random_requests = []

    # read in request files by assignment type
    for proj_request in args.requests:

        request = Request(proj_request, args)
        if request.assignment_type == 'rattaca':
            rattaca_requests.append(request)
        if request.assignment_type == 'hsw_breeders':
            hsw_breeder_requests.append(request)
        if request.assignment_type == 'random':
            random_requests.append(request)

    # save all requests in order of priority:
    # 1. HS West breeders
    # 2. RATTACA projects
    # 3. random assignment - no predictions needed
    all_requests = hsw_breeder_requests + rattaca_requests + random_requests

    # initialize variables to store stats from the eventual best permutation
    best_delta = 0
    best_permutation = None
    best_rfids = None

    # permute the assignment order to all requests
    n_requests = len(all_requests)
    for permuted_order in permutations(range(n_requests)):
        print(f'permuted_order: {permuted_order}')
        proposed_rfids = {}
        permutation_delta = 0
        project_deltas = []

        # for each project request in the current permutation... 
        for current_project in permuted_order:
            print(f'current_project: {current_project}')
            # ...propose animals for assignment to the project 
            # and add the project's delta to the total delta for the permutation
            project_delta, proposal_rfids = all_requests[current_project]\
                .proposal(proposed_rfids) 
            permutation_delta += project_delta
            project_deltas.append(project_delta)
            proposed_rfids[current_project] = proposal_rfids
            print(f'proposal_rfids: {proposal_rfids}')
            print(f'permutation_delta: {permutation_delta}')
            print(f'project_deltas: {project_deltas}')
        # keep the proposed rfids for assignment from the permutation 
        # that maximizes delta
        if permutation_delta >= best_delta:
            best_delta = permutation_delta
            best_permutation = permuted_order
            best_rfids = proposed_rfids
        print(f'current best_permutation: {best_permutation}')
    print(f'final best_permutation: {best_permutation}')
    print(f'best_rfids: {best_rfids}')
    print(f'best_delta: {best_delta}')
    # assign rats to each request in the permuted order that maximizes delta
    for project_index in best_permutation:
        print(project_index)
        request = all_requests[project_index]

        assign_rfids = best_rfids[project_index]
        remove_rfids = [item for sublist in \
                        [v for k,v in best_rfids.items() if k != project_index] \
                            for item in sublist]
        print(assign_rfids)
        print(remove_rfids)

        request.assign(assign_rfids)
        request.remove(remove_rfids)
        

    # for i, request in enumerate(all_requests):
    #     print(f'request: {request.project}')
    #     for project_index in best_permutation:
    #         if i == project_index:
    #             request.assign(best_rfids[i])
    #             print(f'assign: {best_rfids}')
    #         else:
    #             request.remove(best_rfids[i])
    #             print(f'remove: {best_rfids}')


### classes ###

class Request:
    
    def __init__(self, request_file, args):
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadat = json.load(rf)
            self.assignment_type = req_metadat['assignment_type']
            self.project = req_metadat['project']
            self.trait = req_metadat['trait']
            self.n_requested_total = req_metadat['n_rats']['total']
            self.n_requested_males = req_metadat['n_rats']['male']['total']
            self.n_requested_females = req_metadat['n_rats']['female']['total']
            self.n_requested_males_high = req_metadat['n_rats']['male']['high']
            self.n_requested_males_low = req_metadat['n_rats']['male']['low']
            self.n_requested_females_high = req_metadat['n_rats']['female']['high']
            self.n_requested_females_low = req_metadat['n_rats']['female']['low']
            self.n_requested_high = self.n_requested_males_high + self.n_requested_females_high
            self.n_requested_low = self.n_requested_males_low + self.n_requested_females_low

        # initialize dictionaries to hold assigned rats
        self.assigned_males_high = {}
        self.assigned_males_low = {}
        self.assigned_females_high = {}
        self.assigned_females_low = {}
        self.assigned_high = \
            self.assigned_males_high | self.assigned_females_high
        self.assigned_low = \
            self.assigned_males_low | self.assigned_females_low
        self.assigned_males = \
            self.assigned_males_high | self.assigned_males_low
        self.assigned_females = \
        self.assigned_females_high | self.assigned_females_low
        self.assigned_total = self.assigned_males | self.assigned_females
        
        # counters to track the number of assigned rats
        self.n_assigned_males_high = len(self.assigned_males_high)
        self.n_assigned_males_low = len(self.assigned_males_low)
        self.n_assigned_females_high = len(self.assigned_females_high)
        self.n_assigned_females_low = len(self.assigned_females_low)
        self.n_assigned_high = len(self.assigned_high)
        self.n_assigned_low = len(self.assigned_low)
        self.n_assigned_males = len(self.assigned_males)
        self.n_assigned_females = len(self.assigned_females)
        self.n_assigned_total = len(self.assigned_total)

        # counters to track the number of assignments remaining
        self.n_remaining_males_high = \
            self.n_requested_males_high - self.n_assigned_males_high
        self.n_remaining_males_low = \
            self.n_requested_males_low - self.n_assigned_males_low
        self.n_remaining_females_high = \
            self.n_requested_females_high - self.n_assigned_females_high
        self.n_remaining_females_low = \
            self.n_requested_females_low - self.n_assigned_females_low
        self.n_remaining_high = self.n_requested_high - self.n_assigned_high
        self.n_remaining_low = self.n_requested_low - self.n_assigned_low
        self.n_remaining_males = \
            self.n_requested_males - self.n_assigned_males
        self.n_remaining_females = \
            self.n_requested_females - self.n_assigned_females
        self.n_remaining_total = self.n_requested_total - self.n_assigned_total
        
        # initialize the delta for the request
        self.delta = 0

        # read in colony data, predictions, and exclusion list
        colony_df = pd.read_csv(args.colony_dataframe[0], 
                                dtype = {'rfid': str, 'accessid': int})
        predictions_df = pd.read_csv(args.predictions[0], dtype = {'rfid': str})
        
        if args.exclude:
            exclude_rfids = pd.read_csv(args.exclude[0], header = None)\
                .astype(str)[0].values
        else:
            exclude_rfids = []

        # check for duplicate RFIDs in the colony dataframe
        if sum(colony_df.duplicated(subset = ['rfid'])) > 0:
            print('Error: Colony metadata has duplicate RFIDs')
            print('Check the following samples before proceeding:')
            dups = colony_df.duplicated(subset = ['rfid'])
            dups = colony_df[dups]
            for rfid in dups['rfid']:
                print(rfid)
            exit()

        # merge predictions with colony metadata, drop RFIDs to exclude
        rats_metadata = colony_df.merge(predictions_df, on='rfid', how='right')

        # sort by the requested trait, remove NAs
        self.trait_metadata = rats_metadata.sort_values(by = self.trait, 
            axis = 0, ascending = False)\
                .dropna(subset = self.trait, ignore_index=True)

        # add a column of high/low group assignments
        self.trait_metadata[f'{self.trait}_2group'] = \
            self.trait_groups(self.trait, n_groups=2)

        # list of rats available for assignment, ordered by trait prediction
        self.available_rfids = self.trait_metadata['rfid'].tolist()
        self.available_rfids = [rfid for rfid in self.available_rfids 
                                    if rfid not in exclude_rfids]

        # convert df to dictionary with key:RFID and value:(sex, prediction)
        self.trait_metadata = self.trait_metadata[self.trait_metadata['rfid']\
                                                  .isin(self.available_rfids)]
        self.rat_data = self.trait_metadata.\
            set_index('rfid')[['sex', self.trait]].to_dict(orient='index')
        
        # convert the nested dictionaries to tuples
        self.available_rats = \
            {k: (v['sex'], v[self.trait]) for k, v in self.rat_data.items()}
        self.available_males = \
            {k: v for k, v in self.available_rats.items() if v[0] == 'M'}
        self.available_females = \
            {k: v for k, v in self.available_rats.items() if v[0] == 'F'}

    # function to propose an assignment of RFIDs to a project
    #### agnostic to sex, for now
    #### TO DO: keep track of sex: do min/max by sex
    # unavail_rats = either empty, or a list of RFIDs,
    # or a dict of (sex, pred) w/ RFIDs as keys
    def proposal(self, unavail_rats = None):
        
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

        # calculate delta for this iteration: the difference in trait values
        #### TO DO: calculate sex-specific deltas
        min_rat_sex, min_rat_pred = self.available_rats[min_rat]
        max_rat_sex, max_rat_pred = self.available_rats[max_rat]
        print(f'proposal: min_rat: {min_rat}, pred: {min_rat_pred}')
        print(f'proposal: max_rat: {max_rat}, pred: {max_rat_pred}')
        delta_iter = max_rat_pred - min_rat_pred
        print(f'proposal: delta_iter: {delta_iter}')
        # calculate the proposed total delta: the sum of all iterations' 
        # deltas if this proposed iteration were to be added
        proposed_delta = self.delta + delta_iter

        # return proposed delta, selected rats
        return (proposed_delta, [min_rat, max_rat])

    # function to remove assigned rats from the list of available rats
    ## TO DO: remove from M/F, H/L sublists
    def remove(self, rfids_to_remove):
        
        for rfid in rfids_to_remove:
            
            if rfid in self.available_rats:
                self.available_rfids.remove(rfid)
                del self.available_rats[rfid]
                self.n_remaining_total -= 1

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

    # function to assign rats to RATTACA projects based on their genetic predictions
    def assign_rattaca(self, rats_to_assign): # rfids_to_assign = best_rfids

        if isinstance(rats_to_assign, dict):
            max_rat = max(rats_to_assign, key=rats_to_assign.get)
            min_rat = min(rats_to_assign, key=rats_to_assign.get)
        elif isinstance(rats_to_assign, list):
            # get trait values from metadata
            rfids_to_assign = rats_to_assign
            rats_to_assign = {}
            for rfid in rfids_to_assign:
                rats_to_assign[rfid] = self.available_rats[rfid]
            max_rat = max(rats_to_assign, key=rats_to_assign.get)
            min_rat = min(rats_to_assign, key=rats_to_assign.get)
        else:
            raise TypeError('rats_to_assign must be either a list or a dictionary')

        # if the project still needs more assignments, 
        # add the rfids with max and min trait values 
        # to high and low groups, respectively
        # while self.is_satisfied == False:
        n_remaining_males_high = \
            self.n_requested_males_high - self.n_assigned_males_high
        n_remaining_males_low = \
            self.n_requested_males_low - self.n_assigned_males_low
        n_remaining_females_high = \
            self.n_requested_females_high - self.n_assigned_females_high
        n_remaining_females_low = \
            self.n_requested_females_low - self.n_assigned_females_low

            # # get the delta betwen assigned rats
            # min_rat_sex, min_rat_pred = rats_to_assign[min_rat]
            # max_rat_sex, max_rat_pred = rats_to_assign[max_rat]
            # assignment_delta = max_rat_pred - min_rat_pred

            # #### temp for now - remove once sex-specific assignment is adopted
            # self.assigned_low[min_rat] = self.available_rats[min_rat]
            # self.assigned_high[max_rat] = self.available_rats[max_rat]

            # #### unused for now - use this once sex-specific assignment is adopted
            # if min_rat_sex == 'M':
            #     self.assigned_males_low[min_rat] = self.available_rats[min_rat]
            # else:
            #     self.assigned_females_low[min_rat] = self.available_rats[min_rat]
            # if max_rat_sex == 'M':
            #     self.assigned_males_high[max_rat] = self.available_rats[max_rat]
            # else:
            #     self.assigned_females_high[max_rat] = self.available_rats[max_rat]

            # # update the overall delta for the request: add the assignment delta
            # print(f'delta orig: {self.delta}')
            # print(f'assignment_delta: {assignment_delta}')
            # self.delta += assignment_delta
            # print(f'delta new: {self.delta}')

            # # remove assigned rats from further availability
            # self.remove(rfids_to_assign)

    # generic function to check if the request has been fulfilled
    def is_satisfied(self):
        
        if self.assignment_type == 'rattaca':
            return self.is_satisfied_rattaca()
        elif self.assignment_type == 'random':
            return self.is_satisfied_random()
        elif self.assignment_type == 'hsw_breeders':
            return self.is_satisfied_hsw_breeders()
        else:
            raise NameError('Assignment type must be either "rattaca", \
                            "hsw_breeders", or "random". Please check the .json \
                            file for this project request')
            
    # function to check if a prediction-based RATTACA project has 
    # successfully assigned all requested rats
    def is_satisfied_rattaca(self):

        return (
            self.n_assigned_males_high == self.n_requested_males_high and 
            self.n_assigned_males_low == self.n_requested_males_low and
            self.n_assigned_females_high == self.n_requested_females_high and
            self.n_assigned_females_low == self.n_requested_females_low
        )
    
    def is_satisfied_random(self):
        pass

    def is_satisfied_hsw_breeers(self):
        pass


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

    # @property
    # def delta(self):
    #     # raise NotImplementedError
    #     # return self.high - self.low
    #     return self._delta
    
    # @delta.setter
    # def delta(self, value):
    #     self._delta = value


### execute ###

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)