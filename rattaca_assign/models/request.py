import json
import pandas as pd
from datetime import datetime, timedelta
from rattaca_assign.core.base_utils import prep_colony_df

### Request parent class
class Request:
    
    def __init__(self, request_file, args):

        # read in the request metadata from json
        with open(request_file, 'r') as rf:

            req_metadata = json.load(rf)
            self.assignment_type = req_metadata['assignment_type']
            self.project = req_metadata['project']
            self.request_name = req_metadata['request_name']
            self.assignment_name = req_metadata['assignment_name']
            self.trait = req_metadata['trait']
            self.n_requested_total = req_metadata['n_rats']['total']
            self.n_requested_males = req_metadata['n_rats']['male']['total']
            self.n_requested_females = req_metadata['n_rats']['female']['total']
            self.min_age = req_metadata['min_age']
            self.max_age = req_metadata['max_age']
            date_str = str(req_metadata['receive_date']) if req_metadata['receive_date'] is not None else None
            self.receive_date = datetime.strptime(date_str, '%Y%m%d').date() if date_str is not None else None
                
        # read in colony data 
        self.colony_df = prep_colony_df(args)

        # add 'trait' info to dataframe
        self.trait_metadata = prep_colony_df(args)
        self.trait_metadata[[self.trait]] = 0

        # create a list of rats available for assignment, ordered by trait prediction
        self.available_rfids = self.trait_metadata['rfid'].tolist()

        # if only one sex is being requested, remove the opposite sex from availability
        all_males = set(self.trait_metadata[self.trait_metadata['sex'] == 'M']['rfid'].tolist())
        all_females = set(self.trait_metadata[self.trait_metadata['sex'] == 'F']['rfid'].tolist())
        if self.n_requested_males == 0:
            self.available_rfids = [rfid for rfid in self.available_rfids if rfid not in all_males]
        if self.n_requested_females == 0:
            self.available_rfids = [rfid for rfid in self.available_rfids if rfid not in all_females]

        # remove rats from availability that do not meet age restrictions
        if self.receive_date is not None:
            rats_to_remove = []
            for rat in self.available_rfids:
                rat_dob = self.trait_metadata[self.trait_metadata['rfid'] == rat]['dob'].iloc[0]
                rat_dob = datetime.strptime(rat_dob, '%Y-%m-%d').date()
                t_delta = self.receive_date - rat_dob
                shipment_age = t_delta.days
                if shipment_age < self.min_age or shipment_age > self.max_age:
                    rats_to_remove.append(rat)
            for rat in rats_to_remove:
                self.available_rfids.remove(rat)

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
    def assign(self, rats_to_assign=None, remaining_requests=None, fill_request=False):
        '''
        Algorithmically assign RFIDs to a request.

        This is a generic function. Depending on the request type, it will call 
        a type-specific function to assign RFIDs.
        
        args:
            rats_to_assign (optional): a list or dictionary of RFIDs to assign 
                to the request
            remaining_requests (optional): a list of request objects whose 
                assignment has not yet been completed
            assign_remainder (optional): a boolean indicating whether to 
                immediately complete all assignments to the request (without 
                permuting assignments to other requests)
        '''

        if self.assignment_type == 'rattaca':
            return self.assign_rattaca(rats_to_assign)
        elif self.assignment_type == 'random':
            return self.assign_random(rats_to_assign)
        elif self.assignment_type == 'hsw_breeders':
            return self.assign_hsw_breeders(non_breeder_requests=remaining_requests, 
                assign_remainder=fill_request)
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
            return self.is_satisfied_random(sex)
        elif self.assignment_type == 'hsw_breeders':
            return self.is_satisfied_hsw_breeders(sex, family)
        else:
            raise NameError('''Assignment type must be either "rattaca", 
                            "hsw_breeders", or "random". Please check the .json 
                            file for this project request''')
