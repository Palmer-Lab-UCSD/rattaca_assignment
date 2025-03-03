import json
import pandas as pd
from rattaca_assign.utils import prep_colony_df

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
        '''
        Algorithmically assign RFIDs to a request.

        This is a generic function. Depending on the request type, it will call 
        a type-specific function to assign RFIDs.
        
        args:
            rats_to_assign (optional): a list or dictionary of RFIDs to assign 
                to the request
            remaining_requests (optional): a list of request objects whose 
                assignment has not yet been completed
            fill_request (optional): a boolean indicating whether to immediately 
                complete all assignments to the request (without permuting 
                assignments to other requests)
        '''

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
        '''
        Manually assign RFIDs to a request.

        This is a generic function. Depending on the request type, it will call 
        a type-specific function to assign RFIDs.
        
        args:
            rats_to_assign: a list or dictionary of RFIDs to assign 
                to the request
            remaining_requests (optional): a list of request objects whose 
                assignment has not yet been completed
            override (optional): For hsw_breeders requests only. A boolean 
                indicating whether to override constraints set for the breeders
                request. For example, if one male from litter X has already been
                assigned but another is still needed (say to compensate for 
                another litter that produced no males), setting override = False
                will result in an error due to request constraints, but 
                override = True will successfully assign all IDs provided by 
                rats_to_assign, regardless of constraints. Using override = False
                is a good strategy for assigning breeders ad-hoc to replace 
                original assignments that may be dead or missing, while ensuring
                desired constraints are still observed.
        '''

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
        '''
        Test whether a request has been fulfilled (all assignments to the 
        request have been satisfied).

        Any combination of sexes, group classifications, or family IDs can be 
        used. Whichever parameters are set, the funciton will test whether the
        total number of rats meeting those criteria have been successfully 
        assigned to the request. 

        For example, Request.is_satisfied() tests if all assignments (to all
        groups) have been completed for the request. 
        Request.is_satisfied(sex = 'M', group = 'low') tests if all assignments 
        for high-group males have been completed. 

        This is a generic function. Depending on the request type, it will call 
        a type-specific function to test request fulfillment.
        
        Args:
            sex (optional): Tests if the total number of rats of a given 
                sex have been assigned
            group (optional): Tests if the total number of rats from a given 
                group classification (high vs. low) have been assigned
            family (optional): Tests if the total number of rats allowed from a 
                given breeder pair have been assigned

        Returns:
            A boolean indicating request fulfillment under the constraints set 
            by the input parameters.
        '''
        
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
