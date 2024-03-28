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
    n_requests = len(args.requests)
    for permuted_order in permutations(range(n_requests)):
        
        tmp_unavailable_rfids = []
        permutation_delta = 0

        # for each project request in the current permutation 
        for current_project in permuted_order:
            # get the project's overall delta & set of animals assigned 
            # to the project
            project_delta, proposed_rfids = all_requests[current_project]\
                .proposal(tmp_unavailable_rfids) ### edit proposal() ###
            permutation_delta += project_delta
            tmp_unavailable_rfids.append(proposed_rfids) ### edit append() ###

        if permutation_delta > best_delta:
            best_delta = permutation_delta
            best_permutation = permuted_order
            best_rfids = tmp_unavailable_rfids
    
    # assign rats to each request in the permuted order that maximizes delta
    for i, request in enumerate(all_requests):

        for project_index in best_permutation:
            if i == project_index:
                request.update(best_rfids[i])
            else:
                request.remove(best_rfids[i])


### classes ###

class Request:
    
    def __init__(self, request_file, args):
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadat = json.load(rf)
            self.project = req_metadat['project']
            self.trait = req_metadat['trait']
            self.assignment_type = req_metadat['assignment_type']
            self.total_rats = req_metadat['n_rats']['total']
            self.total_males = req_metadat['n_rats']['male']['total']
            self.total_females = req_metadat['n_rats']['female']['total']
            self.n_males_high = req_metadat['n_rats']['male']['high']
            self.n_males_low = req_metadat['n_rats']['male']['low']
            self.n_females_high = req_metadat['n_rats']['female']['high']
            self.n_females_low = req_metadat['n_rats']['female']['low']

        # read in colony data, predictions, and exclusion list
        colony_df = pd.read_csv(args.colony_dataframe[0], 
                                dtype = {'rfid': str, 'accessid': int})
        predictions_df = pd.read_csv(args.predictions[0], dtype = {'rfid': str})
        
        if args.exclude:
            exclude_rfids = pd.read_csv(args.exclude[0], header = None)\
                .astype(str)[0].values
        
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
        if args.exclude: 
            rats_metadata = rats_metadata[~rats_metadata['rfid']\
                .isin(exclude_rfids)]

        # sort by the requested trait, remove NAs
        self.rats_metadata = rats_metadata.sort_values(by = self.trait, 
            axis = 0, ascending = False)\
                .dropna(subset = self.trait, ignore_index=True)
        # print(f'{ len(self.rats_metadata)} rats are available for selection')
        
        # list of trait predictions from which to pull rats
        self.available_rfids = self.rats_metadata['rfid'].tolist()
        
        self.available_trait_vals = dict(zip(self.rats_metadata['rfid'], 
                                             self.rats_metadata[self.trait]))
        self.available_sex_vals = dict(zip(self.rats_metadata['rfid'], 
                                             self.rats_metadata['sex']))

    # function to propose an assignment of RFIDs to a project
    def proposal(self, unavailable_rfids):
        
        # initial indexes to start search for extreme available rats
        i = 0    # first element is the rat with the max trait value
        j = len(self.available_rfids) - 1 # last element has the min trait value

        # select the most extreme high & low available rats
        while self.available_rfids[i] in unavailable_rfids:
            i += 1
        while self.available_rfids[j] in unavailable_rfids:
            j -= 1
        max_rat = self.available_rfids[i]
        min_rat = self.available_rfids[j]

        # calculate delta for this iteration: the difference in trait values
        delta_iter = self.available_trait_vals[max_rat] - \
            self.available_trait_vals[min_rat]

        # calculate the proposed total delta: the sum of all iterations' deltas,
        # if this proposed iteration were to be added
        proposed_delta = delta_iter + self.delta

        # return proposed delta, selected rats
        return (proposed_delta, [max_rat, min_rat])

    # function to remove assigned rats from the list of available rats
    def remove(self, rfids_to_remove):
        
        for rfid in rfids_to_remove:
            self.available_rfids.remove(rfid)
            del self.available_trait_vals[rfid]
            del self.available_sex_vals[rfid]

    def update(self, rfids_to_assign): # rfids_to_assign = best_rfids
        raise NotImplementedError    

    def is_satisfied(self):
        raise NotImplementedError

    @property
    def delta(self):
        # raise NotImplementedError
        # return self.high - self.low
        return 0


### execute ###

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)