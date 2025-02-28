import json
import pandas as pd
from rattaca_assignment.utils import prep_colony_df

### Request parent class
class Request:
    
    def __init__(self, request_file, args):

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


    # function to remove assigned rats from the list of available rats
    def remove(self, rats_to_remove, remaining_requests = None):
        '''
        Remove assigned rats from the list of available rats

        Args:
            rats_to_remove: A list or dictionary of RFIDs to remove
            remaining_requests: A list of requests that still have outstanding 
                assignments to be filled
        '''
        
        if isinstance(rats_to_remove, list):
            rfids_to_remove = rats_to_remove
        elif isinstance(rats_to_remove, dict):
            rfids_to_remove = list(rats_to_remove.keys())
        else:
            raise TypeError('rats_to_remove must be either a list or a dictionary')

        for rfid in rfids_to_remove:
            
            if rfid in self.available_rfids:
                self.available_rfids.remove(rfid) # python list.remove() method, not RATTACA remove() function 
        
        if self.assignment_type == 'random':
            self.update_randproj(remaining_requests = remaining_requests)
        
        if self.assignment_type == 'hsw_breeders':
            self.update_breeders(non_breeder_requests = remaining_requests)

    # generic function to assign rats to a project, remove them from availability
    def assign(self, rats_to_assign=None, remaining_requests=None, fill_request=None):
        '''Algorithmically assign RFIDs to a request'''

        if self.assignment_type == 'rattaca':
            return self.assign_rattaca(rats_to_assign)
        elif self.assignment_type == 'random':
            return self.assign_random(rats_to_assign)
        elif self.assignment_type == 'hsw_breeders':
            return self.assign_hsw_breeders(non_breeder_requests=remaining_requests, fill_request=fill_request)
        else:
            raise NameError('Assignment type must be either "rattaca", \
                            "hsw_breeders", or "random". Please check the .json \
                            file for this project request')

    def assign_manual(self, rats_to_assign, override=None):
        '''Manually assign RFIDs to a request'''

        if self.assignment_type == 'rattaca':
            return self.assign_manual_rattaca(rats_to_assign)
        elif self.assignment_type == 'random':
            return self.assign_manual_random(rats_to_assign)
        elif self.assignment_type == 'hsw_breeders':
            return self.assign_manual_hsw_breeders(rats_to_assign, override)
        else:
            raise NameError('Assignment type must be either "rattaca", \
                            "hsw_breeders", or "random". Please check the .json \
                            file for this project request')

    # generic function to check if the request has been fulfilled
    def is_satisfied(self, sex=None, group=None, family=None):
        
        if self.assignment_type == 'rattaca':
            return self.is_satisfied_rattaca(sex, group)
        elif self.assignment_type == 'random':
            return self.is_satisfied_random(sex, family)
        elif self.assignment_type == 'hsw_breeders':
            return self.is_satisfied_hsw_breeders(sex, family)
        else:
            raise NameError('''Assignment type must be either "rattaca", 
                            "hsw_breeders", or "random". Please check the .json 
                            file for this project request''')
