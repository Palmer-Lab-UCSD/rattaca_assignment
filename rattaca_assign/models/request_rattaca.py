'''
RATTACA request class for assignment based on genetic predictions.
'''

import json
import numpy as np
import pandas as pd
from rattaca.models.request import Request


class RATTACA(Request):
    '''
    Request subclass for RATTACA projects.
    
    This class handles assignment of rats to projects based on their 
    genetic predictions for specific traits.
    '''
    
    def __init__(self, request_file, args):
        '''
        Initialize a RATTACA request.
        
        Args:
            request_file: Path to JSON request file with request_type='rattaca'
            args: Command line arguments.
        '''
        print('Initializing RATTACA request')

        # inherit attributes from the parent Request class
        super().__init__(request_file, args)
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:
            req_metadata = json.load(rf)
            self.n_requested_males_high = req_metadata['n_rats']['male']['high']
            self.n_requested_males_low = req_metadata['n_rats']['male']['low']
            self.n_requested_females_high = req_metadata['n_rats']['female']['high']
            self.n_requested_females_low = req_metadata['n_rats']['female']['low']
            self.n_requested_high = self.n_requested_males_high + self.n_requested_females_high
            self.n_requested_low = self.n_requested_males_low + self.n_requested_females_low

        # initialize dictionaries to hold assigned rats
        self.assigned_males_high = {}
        self.assigned_males_low = {}
        self.assigned_females_high = {}
        self.assigned_females_low = {}
                
        # initialize the delta for the request
        self.delta = 0

        # merge trait predictions into colony data
        self.setup_trait_metadata(args)

        print('Finished initializing RATTACA request')


    # function to manually assign rats to a projects as needed
    def assign_manual_rattaca(self, rfids_to_assign):
        '''
        Manually assign RFIDs to a RATTACA request.
        
        Args:
            rfids_to_assign: A list of desired RFIDs to be assigned to the project            
        '''
        if not isinstance(rfids_to_assign, list):
            raise TypeError('rfids_to_assign must be a list')
        
        for rfid in rfids_to_assign:
            if rfid not in self.available_rfids:
                print(f'RFID {rfid} is not available for assignment, skipping')
                continue
                
            # Get rat information
            rat_data = self.trait_metadata[self.trait_metadata['rfid'] == rfid].iloc[0]
            rat_sex = rat_data['sex']
            rat_group = rat_data[f'{self.trait}_2group']
            rat_pred = rat_data[self.trait]
            
            # determine assignment destination
            if rat_sex == 'M' and rat_group == 'high':
                if not self.is_satisfied_rattaca('M', 'high'):
                    self.assigned_males_high[rfid] = (rat_sex, rat_pred, rat_group)
                    self.remove([rfid])
                else:
                    print(f'Cannot assign {rfid}: Male/high group already satisfied')
            elif rat_sex == 'M' and rat_group == 'low':
                if not self.is_satisfied_rattaca('M', 'low'):
                    self.assigned_males_low[rfid] = (rat_sex, rat_pred, rat_group)
                    self.remove([rfid])
                else:
                    print(f'Cannot assign {rfid}: Male/low group already satisfied')
            elif rat_sex == 'F' and rat_group == 'high':
                if not self.is_satisfied_rattaca('F', 'high'):
                    self.assigned_females_high[rfid] = (rat_sex, rat_pred, rat_group)
                    self.remove([rfid])
                else:
                    print(f'Cannot assign {rfid}: Female/high group already satisfied')
            elif rat_sex == 'F' and rat_group == 'low':
                if not self.is_satisfied_rattaca('F', 'low'):
                    self.assigned_females_low[rfid] = (rat_sex, rat_pred, rat_group)
                    self.remove([rfid])
                else:
                    print(f'Cannot assign {rfid}: Female/low group already satisfied')
        
        # update available rats list
        self._update_available_rats()
        

    # function to classify RFIDs into high/low groups
    def trait_groups(self, trait, n_groups=2):
        '''
        Classify group values based on trait predictions. This does not
        'officially' assign rats to a group, rather classifies whether a rat 
        belongs into high/low groups (or another quantile) based on its 
        predicted trait value. Samples are drawn from these groups using
        assign_rattaca().
        
        Args:
            trait: The trait to group by.
            n_groups: The number of groups to create.
            
        Returns:
            A list of group classifications.
        '''
        df = self.trait_metadata
        preds = df[trait].tolist()
        trait_quantiles = np.quantile(a=preds, q=np.linspace(0, 1, n_groups+1))
        
        if n_groups == 2:
            def get_2group(val, quantiles):
                if val >= quantiles[1]:
                    return 'high'
                else:
                    return 'low'
            groups = [get_2group(pred, trait_quantiles) for pred in preds]
        
        elif n_groups == 3:
            def get_3group(val, quantiles):
                if val < quantiles[1]:
                    return 'low'
                elif val >= quantiles[2]:
                    return 'high'
                else:
                    return 'mid'
            groups = [get_3group(pred, trait_quantiles) for pred in preds]
        
        else:
            groups = np.digitize(preds, trait_quantiles)
            
        return groups


    @property
    def available_rats(self):
        '''
        A dictionary of all rats currently available for assignment.
        
        Returns:
            A dictionary mapping RFIDs to (sex, prediction, group) tuples.
        '''
        # update trait metadata to include only currently available rats
        current_metadata = self.trait_metadata[
            self.trait_metadata['rfid'].isin(self.available_rfids)
        ]
        
        # convert the metadata df to a dictionary with key:RFID and 
        # value:(sex, prediction, group)
        current_metadata = current_metadata.set_index('rfid')[
            ['sex', self.trait, f'{self.trait}_2group']
        ].to_dict(orient='index')
        
        available_rats = {
            k: (v['sex'], v[self.trait], v[f'{self.trait}_2group']) 
            for k, v in current_metadata.items()
        }
        
        return available_rats
    
    @property
    def available_males(self):
        '''A dictionary of all male rats available for assignment.'''
        return {k: v for k, v in self.available_rats.items() if v[0] == 'M'}
    
    @property
    def available_females(self):
        '''A dictionary of all female rats available for assignment.'''
        return {k: v for k, v in self.available_rats.items() if v[0] == 'F'}

    @property
    def assigned_rats(self):
        '''
        A dictionary of all rats assigned to the project.
        
        Returns:
            A dictionary with different groupings of assigned rats.
        '''
        assigned_high = {**self.assigned_males_high, **self.assigned_females_high}
        assigned_low = {**self.assigned_males_low, **self.assigned_females_low}
        assigned_males = {**self.assigned_males_high, **self.assigned_males_low}
        assigned_females = {**self.assigned_females_high, **self.assigned_females_low}
        assigned_total = {**assigned_males, **assigned_females}
        
        return {
            'assigned_high': assigned_high, 
            'assigned_low': assigned_low,
            'assigned_males': assigned_males, 
            'assigned_females': assigned_females,
            'assigned_total': assigned_total
        }

    @property
    def n_assigned(self):
        '''
        Get counts of assigned rats by category.
        
        Returns:
            A dictionary with counts of assigned rats by category.
        '''
        n_assigned_males_high = len(self.assigned_males_high)
        n_assigned_males_low = len(self.assigned_males_low)
        n_assigned_females_high = len(self.assigned_females_high)
        n_assigned_females_low = len(self.assigned_females_low)
        n_assigned_high = len(self.assigned_rats['assigned_high'])
        n_assigned_low = len(self.assigned_rats['assigned_low'])
        n_assigned_males = len(self.assigned_rats['assigned_males'])
        n_assigned_females = len(self.assigned_rats['assigned_females'])
        n_assigned_total = len(self.assigned_rats['assigned_total'])

        return {
            'n_assigned_males_high': n_assigned_males_high,
            'n_assigned_males_low': n_assigned_males_low,
            'n_assigned_females_high': n_assigned_females_high,
            'n_assigned_females_low': n_assigned_females_low,
            'n_assigned_high': n_assigned_high,
            'n_assigned_low': n_assigned_low,
            'n_assigned_males': n_assigned_males,
            'n_assigned_females': n_assigned_females,
            'n_assigned_total': n_assigned_total
        }

    @property
    def n_remaining(self):
        '''
        Get counts of remaining assignments by category.
        
        Returns:
            A dictionary with counts of remaining assignments by category.
        '''
        n_remaining_males_high = max(0, self.n_requested_males_high - self.n_assigned['n_assigned_males_high'])
        n_remaining_males_low = max(0, self.n_requested_males_low - self.n_assigned['n_assigned_males_low'])
        n_remaining_females_high = max(0, self.n_requested_females_high - self.n_assigned['n_assigned_females_high'])
        n_remaining_females_low = max(0, self.n_requested_females_low - self.n_assigned['n_assigned_females_low'])
        n_remaining_high = max(0, self.n_requested_high - self.n_assigned['n_assigned_high'])
        n_remaining_low = max(0, self.n_requested_low - self.n_assigned['n_assigned_low'])
        n_remaining_males = max(0, self.n_requested_males - self.n_assigned['n_assigned_males'])
        n_remaining_females = max(0, self.n_requested_females - self.n_assigned['n_assigned_females'])
        n_remaining_total = max(0, self.n_requested_total - self.n_assigned['n_assigned_total'])

        return {
            'n_remaining_males_high': n_remaining_males_high,
            'n_remaining_males_low': n_remaining_males_low,
            'n_remaining_females_high': n_remaining_females_high,
            'n_remaining_females_low': n_remaining_females_low,
            'n_remaining_high': n_remaining_high,
            'n_remaining_low': n_remaining_low,
            'n_remaining_males': n_remaining_males,
            'n_remaining_females': n_remaining_females,
            'n_remaining_total': n_remaining_total
        }
        

    def is_satisfied_rattaca(self, sex=None, group=None):
        '''
        Check if a RATTACA request is satisfied.
        
        Args:
            sex: Optional sex filter (M/F).
            group: Optional group filter (high/low).
            
        Returns:
            Boolean indicating if request is satisfied.
        '''
        n_remaining = self._get_n_remaining(sex, group)
        return n_remaining == 0
        

    # function to count how many assignments remain per group
    def _get_n_remaining(self, sex=None, group=None):
        '''Get the number of assignments that remain to be filled 
        for a specific group.'''
        if sex is None and group is None:
            return self.n_remaining['n_remaining_total']
        if sex == 'M' and group is None:
            return self.n_remaining['n_remaining_males']
        elif sex == 'F' and group is None:
            return self.n_remaining['n_remaining_females']
        if group == 'high' and sex is None:
            return self.n_remaining['n_remaining_high']
        if group == 'low' and sex is None:
            return self.n_remaining['n_remaining_low']
        if sex == 'M' and group == 'high':
            return self.n_remaining['n_remaining_males_high']
        elif sex == 'M' and group == 'low':
            return self.n_remaining['n_remaining_males_low']
        elif sex == 'F' and group == 'high':
            return self.n_remaining['n_remaining_females_high']
        elif sex == 'F' and group == 'low':
            return self.n_remaining['n_remaining_females_low']


    def setup_trait_metadata(self, args):
        '''
        Set up trait metadata for a RATTACA request.
        
        Args:
            args: Command line arguments.
        '''

        if args.exclude:
            exclude_rfids = pd.read_csv(args.exclude[0], header=None).astype(str)[0].values
        else:
            exclude_rfids = []

        # read in predictions
        predictions_df = pd.read_csv(args.predictions[0], dtype={'rfid': str})

        # merge predictions with colony metadata
        all_metadata = self.colony_df.merge(
            predictions_df[['rfid', self.trait]], on='rfid', how='right')

        # sort by the requested trait, remove NAs
        self.trait_metadata = all_metadata.sort_values(
            by=self.trait, axis=0, ascending=False
        ).dropna(subset=self.trait, ignore_index=True)

        # add a column of high/low group assignments
        self.trait_metadata[f'{self.trait}_2group'] = self.trait_groups(self.trait, n_groups=2)

        # drop rats from the exclude list
        self.trait_metadata = self.trait_metadata[~self.trait_metadata['rfid'].isin(exclude_rfids)]

        # drop any groups that aren't requested
        if self.n_requested_males == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata['sex'] != 'M']
        if self.n_requested_females == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata['sex'] != 'F']
        if self.n_requested_high == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata[f'{self.trait}_2group'] != 'high']
        if self.n_requested_low == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata[f'{self.trait}_2group'] != 'low']

        # create a list of rats available for assignment, ordered by trait prediction
        self.available_rfids = self.trait_metadata['rfid'].tolist()


    def proposal(self, unavail_rats=None):
        '''
        Propose rats for assignment to this RATTACA request.
        
        Args:
            unavail_rats: Optional dictionary of unavailable RFIDs.
            
        Returns:
            Tuple of (delta, proposed_rfids).
        '''
        if unavail_rats is None:
            unavail_rfids = []
        elif isinstance(unavail_rats, list):
            unavail_rfids = unavail_rats
        elif isinstance(unavail_rats, dict):
            unavail_rfids = []
            for rfids_list in unavail_rats.values():
                unavail_rfids.extend(rfids_list)
        else:
            raise TypeError('unavail_rats must be either a list, a dictionary, or None')
        
        # initial indexes to start search for extreme available rats
        i = 0    # first element is the rat with the max trait value
        j = len(self.available_rfids) - 1 # last element has the min trait value

        # make sure enough rats are available to make a proposal
        if len(self.available_rfids) < 2 or all(rfid in unavail_rfids for rfid in self.available_rfids):
            return 0, []
            
        # find rats that are currently available
        available = [rfid for rfid in self.available_rfids if rfid not in unavail_rfids]
        if len(available) < 2:
            return 0, []
            
        # get the highest and lowest predictions from available rats
        available_df = self.trait_metadata[self.trait_metadata['rfid'].isin(available)]
        max_rat = available_df.iloc[0]['rfid']  # first row (highest prediction)
        min_rat = available_df.iloc[-1]['rfid']  # last row (lowest prediction)
        
        # get prediction values
        max_rat_pred = available_df[available_df['rfid'] == max_rat][self.trait].iloc[0]
        min_rat_pred = available_df[available_df['rfid'] == min_rat][self.trait].iloc[0]
        
        # calculate delta for this proposal
        delta_iter = max_rat_pred - min_rat_pred
        
        # calculate the proposed total delta if this proposal were to be executed
        proposed_delta = self.delta + delta_iter

        # return proposed delta and selected rats
        return proposed_delta, [min_rat, max_rat]


    def assign_rattaca(self, rfids_to_assign):
        '''
        Assign rats to a RATTACA project.
        
        Args:
            rfids_to_assign: List of RFIDs to assign, ordered [min_rat, max_rat].
        '''

        if not isinstance(rfids_to_assign, list) or len(rfids_to_assign) != 2:
            raise TypeError('rfids_to_assign must be a list with two RFID elements')
        
        # make sure rats are available
        available_rfids = set(self.available_rfids)
        for rfid in rfids_to_assign:
            if rfid not in available_rfids:
                print(f'RFID {rfid} is not available for assignment, skipping')
                
        # get rat info
        min_rat = rfids_to_assign[0]
        max_rat = rfids_to_assign[1]
        
        min_rat_data = self.trait_metadata[self.trait_metadata['rfid'] == min_rat].iloc[0]
        max_rat_data = self.trait_metadata[self.trait_metadata['rfid'] == max_rat].iloc[0]
        
        min_rat_sex = min_rat_data['sex']
        min_rat_group = min_rat_data[f'{self.trait}_2group']
        min_rat_pred = min_rat_data[self.trait]
        
        max_rat_sex = max_rat_data['sex']
        max_rat_group = max_rat_data[f'{self.trait}_2group']
        max_rat_pred = max_rat_data[self.trait]
        
        # assign the high rat
        if max_rat_sex == 'M':
            if not self.is_satisfied_rattaca('M', 'high'):
                self.assigned_males_high[max_rat] = (max_rat_sex, max_rat_pred, max_rat_group)
                self.remove([max_rat])
            else:
                print(f'Cannot assign {max_rat}: Male/high group already satisfied')
        elif max_rat_sex == 'F':
            if not self.is_satisfied_rattaca('F', 'high'):
                self.assigned_females_high[max_rat] = (max_rat_sex, max_rat_pred, max_rat_group)
                self.remove([max_rat])
            else:
                print(f'Cannot assign {max_rat}: Female/high group already satisfied')
        
        # assign the low rat
        if min_rat_sex == 'M':
            if not self.is_satisfied_rattaca('M', 'low'):
                self.assigned_males_low[min_rat] = (min_rat_sex, min_rat_pred, min_rat_group)
                self.remove([min_rat])
            else:
                print(f'Cannot assign {min_rat}: Male/low group already satisfied')
        elif min_rat_sex == 'F':
            if not self.is_satisfied_rattaca('F', 'low'):
                self.assigned_females_low[min_rat] = (min_rat_sex, min_rat_pred, min_rat_group)
                self.remove([min_rat])
            else:
                print(f'Cannot assign {min_rat}: Female/low group already satisfied')
        
        # update delta
        self.delta += max_rat_pred - min_rat_pred
        
        # update available rats list
        self._update_available_rats()
        

    def _update_available_rats(self):
        '''Update available rats list based on satisfied groups.'''
        if self.is_satisfied_rattaca('M', 'high'):
            high_male_rfids = self.trait_metadata[
                (self.trait_metadata['sex'] == 'M') & 
                (self.trait_metadata[f'{self.trait}_2group'] == 'high')
            ]['rfid'].tolist()
            for rfid in high_male_rfids:
                if rfid in self.available_rfids:
                    self.available_rfids.remove(rfid)
                    
        if self.is_satisfied_rattaca('M', 'low'):
            low_male_rfids = self.trait_metadata[
                (self.trait_metadata['sex'] == 'M') & 
                (self.trait_metadata[f'{self.trait}_2group'] == 'low')
            ]['rfid'].tolist()
            for rfid in low_male_rfids:
                if rfid in self.available_rfids:
                    self.available_rfids.remove(rfid)
                    
        if self.is_satisfied_rattaca('F', 'high'):
            high_female_rfids = self.trait_metadata[
                (self.trait_metadata['sex'] == 'F') & 
                (self.trait_metadata[f'{self.trait}_2group'] == 'high')
            ]['rfid'].tolist()
            for rfid in high_female_rfids:
                if rfid in self.available_rfids:
                    self.available_rfids.remove(rfid)
                    
        if self.is_satisfied_rattaca('F', 'low'):
            low_female_rfids = self.trait_metadata[
                (self.trait_metadata['sex'] == 'F') & 
                (self.trait_metadata[f'{self.trait}_2group'] == 'low')
            ]['rfid'].tolist()
            for rfid in low_female_rfids:
                if rfid in self.available_rfids:
                    self.available_rfids.remove(rfid)
