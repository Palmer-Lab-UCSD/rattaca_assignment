'''
RATTACA request class for assignment based on genetic predictions.
'''

import json
import numpy as np
import pandas as pd
from datetime import datetime
from rattaca_assign.models.request import Request


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

        # inherit attributes from the parent Request class
        super().__init__(request_file, args)
        
        # read in the request metadata from json
        with open(request_file, 'r') as rf:
            req_metadata = json.load(rf)
            self.trait = req_metadata['trait']
            self.project = req_metadata['project']
            self.n_requested_males_high = req_metadata['n_rats']['male']['high']
            self.n_requested_males_low = req_metadata['n_rats']['male']['low']
            self.n_requested_females_high = req_metadata['n_rats']['female']['high']
            self.n_requested_females_low = req_metadata['n_rats']['female']['low']
            self.n_requested_high = self.n_requested_males_high + self.n_requested_females_high
            self.n_requested_low = self.n_requested_males_low + self.n_requested_females_low
            self.max_per_sex = req_metadata['max_per_sex']
            self.min_age = req_metadata['min_age']
            self.max_age = req_metadata['max_age']
            date_str = str(req_metadata['receive_date']) if req_metadata['receive_date'] is not None else None
            self.receive_date = datetime.strptime(date_str, "%Y%m%d").date() if date_str is not None else None
        
        # initialize dictionaries to hold assigned rats
        self.assigned_males_high = {}
        self.assigned_males_low = {}
        self.assigned_females_high = {}
        self.assigned_females_low = {}
                
        # initialize the delta for the request
        self.delta = 0

        # merge trait predictions into colony data
        self.setup_trait_metadata(args)

        ## identify families with rats available for assignment by group
        self.all_fams = self.colony_df['breederpair'].unique().tolist()

        md_males_high = self.trait_metadata\
            [(self.trait_metadata['sex'] == 'M') & \
             (self.trait_metadata[f'{self.trait}_group'] == 'high')]
        md_males_low = self.trait_metadata\
            [(self.trait_metadata['sex'] == 'M') & \
             (self.trait_metadata[f'{self.trait}_group'] == 'low')]
        md_females_high = self.trait_metadata\
            [(self.trait_metadata['sex'] == 'F') & \
             (self.trait_metadata[f'{self.trait}_group'] == 'high')]
        md_females_low = self.trait_metadata\
            [(self.trait_metadata['sex'] == 'F') & \
             (self.trait_metadata[f'{self.trait}_group'] == 'low')]
        
        self.fams_with_high_males = md_males_high['breederpair'].unique()
        self.fams_with_low_males = md_males_low['breederpair'].unique()
        self.fams_with_high_females = md_females_high['breederpair'].unique()
        self.fams_with_low_females = md_females_low['breederpair'].unique()

        self.fams_without_males = list(set(self.all_fams) - \
            set(self.fams_with_high_males) - set(self.fams_with_low_males)) 
        self.fams_without_females = list(set(self.all_fams) - \
            set(self.fams_with_high_females) - set(self.fams_with_low_females)) 

        print(f'Initialized RATTACA request: {self.trait} for {self.project}')


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
            rat_group = rat_data[f'{self.trait}_group']
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
            ['sex', self.trait, f'{self.trait}_group']
        ].to_dict(orient='index')
        
        available_rats = {
            k: (v['sex'], v[self.trait], v[f'{self.trait}_group']) 
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
        

    # property to keep track of all breeder pairs that have been assigned to 
    # open assignment groups
    @property
    def assigned_fams(self):
        '''
        Get all breederpairs that have contributed rats to the request.
        
        Returns:
            A dictionary with breederpair IDs of pairs with offspring that have
            been assigned as to the request, grouped by sex and assignment group. 
        '''

        assigned_fams_m_high = {}
        assigned_fams_m_low = {}
        assigned_fams_f_high = {}
        assigned_fams_f_low = {}

        # update trait metadata to include only breederpairs that have assigned rats
        current_metadata = self.trait_metadata[
            self.trait_metadata['rfid'].isin(self.assigned_rats['assigned_total'])]
        
        # subset assigned rats by assignment group
        current_md_high = current_metadata[current_metadata[f'{self.trait}_group'] == 'high']
        current_md_low = current_metadata[current_metadata[f'{self.trait}_group'] == 'low']
        current_md_m_high = current_md_high[current_md_high['sex'] == 'M']
        current_md_m_low = current_md_low[current_md_low['sex'] == 'M']
        current_md_f_high = current_md_high[current_md_high['sex'] == 'F']
        current_md_f_low = current_md_low[current_md_low['sex'] == 'F']
        
        # convert metadata dfs to dictionaries
        # with key:breederpair and value:(rfid, sex)
        if self.n_assigned['n_assigned_males_high'] > 0:
            assigned_fams_m_high = current_md_m_high\
                .groupby('breederpair')[['rfid', 'sex']]\
                .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        if self.n_assigned['n_assigned_males_low'] > 0:
            assigned_fams_m_low = current_md_m_low\
                .groupby('breederpair')[['rfid', 'sex']]\
                .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        if self.n_assigned['n_assigned_females_high'] > 0:
            assigned_fams_f_high = current_md_f_high\
                .groupby('breederpair')[['rfid', 'sex']]\
                .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        if self.n_assigned['n_assigned_females_low'] > 0:
            assigned_fams_f_low = current_md_f_low\
                .groupby('breederpair')[['rfid', 'sex']]\
                .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()

        out = {'M_high' : assigned_fams_m_high,
               'M_low' : assigned_fams_m_low,
               'F_high' : assigned_fams_f_high,
               'F_low' : assigned_fams_f_low}
            
        return out

    # property to keep track of all breeder pairs currently available
    # to sample for assignment
    @property
    def available_fams(self):
        '''
        Get all breederpairs from which rats can be, but have not yet 
        been assigned.
        
        Returns:
            A dictionary with breederpair IDs for pairs with offspring that are 
            still available for assignment to the request.
        '''

        assigned_fams_m_high = list(self.assigned_fams\
            ['M_high'].keys())
        assigned_fams_m_low = list(self.assigned_fams\
            ['M_low'].keys())
        assigned_fams_f_high = list(self.assigned_fams\
            ['F_high'].keys())
        assigned_fams_f_low = list(self.assigned_fams\
            ['F_low'].keys())

        available_fams_m_high = set(self.fams_with_high_males) - \
            set(assigned_fams_m_high) - set(self.fams_without_males)
        available_fams_m_low = set(self.fams_with_low_males) - \
            set(assigned_fams_m_low) - set(self.fams_without_males)
        available_fams_f_high = set(self.fams_with_high_females) - \
            set(assigned_fams_f_high) - set(self.fams_without_females)
        available_fams_f_low = set(self.fams_with_low_females) - \
            set(assigned_fams_f_low) - set(self.fams_without_females)

        current_md_m_high = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_m_high)]
        current_md_m_high = current_md_m_high\
            [(current_md_m_high['sex'] == 'M') & \
             (current_md_m_high[f'{self.trait}_group'] == 'high')]
        current_md_m_low = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_m_low)]
        current_md_m_low = current_md_m_low\
            [(current_md_m_low['sex'] == 'M') & \
             (current_md_m_low['{self.trait}_group'] == 'low')]
        current_md_f_high = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_f_high)]
        current_md_f_high = current_md_f_high\
            [(current_md_f_high['sex'] == 'F') & \
             (current_md_f_high['{self.trait}_group'] == 'high')]
        current_md_f_low = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_f_low)]
        current_md_f_low = current_md_f_low\
            [(current_md_f_low['sex'] == 'F') & \
             (current_md_f_low['{self.trait}_group'] == 'low')]

        available_fams_m_high = \
            current_md_m_high.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        available_fams_m_low = \
            current_md_m_low.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        available_fams_f_high = \
            current_md_f_high.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        available_fams_f_low = \
            current_md_f_low.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()

        out = {'M_high': available_fams_m_high,
               'M_low': available_fams_m_low,
               'F_high': available_fams_f_high,
               'F_low': available_fams_f_low}

        return out


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


    # count how many remaining assignments are allowed per family x sex
    def _get_n_remaining_per_fam(self, sex, fam):
        '''Get the number of assignments that are still allowed to be drawn from
        a specific family for a given sex'''

        # get the number of rats (of the desired sex) already assigned to the family
        n_sibs_assigned = 0
        for key, breederpair_dict in self.assigned_fams.items():
            if sex in key:
                if fam in breederpair_dict:
                    n_sibs_assigned += len(breederpair_dict[fam])

        # available remainder is the max rats allowed per sex per fam minus the number already assigned
        n_remaining = self.max_per_sex - n_sibs_assigned
        
        return n_remaining


    def _rfid_metadata(self, rfid):
        '''Get relevant metadata for a specific RFID'''

        id_df = self.trait_metadata[self.trait_metadata['rfid'] == rfid].iloc[0]
        
        id_sex = id_df['sex']
        id_fam = id_df['breederpair']
        id_group = id_df[f'{self.trait}_group']
        id_pred = id_df[self.trait]

        # get RFIDs of same-sex siblings
        sib_ids = self.trait_metadata\
            [(self.trait_metadata['breederpair'] == id_fam) & 
             (self.trait_metadata['sex'] == id_sex) & 
             (self.trait_metadata['rfid'] != rfid)]['rfid'].tolist()
        out = {'rfid': rfid,
            'sex': id_sex,
            'group': id_group,
            'pred': id_pred,
            'fam': id_fam,
            'sibs': sib_ids}
        
        return out

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
            predictions_df[['rfid', self.trait]], on='rfid', how='inner')

        # sort by the requested trait, remove NAs
        self.trait_metadata = all_metadata.sort_values(
            by=self.trait, axis=0, ascending=False
        ).dropna(subset=self.trait, ignore_index=True)

        # add a column of high/low group assignments
        self.trait_metadata[f'{self.trait}_group'] = self.trait_groups(self.trait, n_groups=2)

        # drop rats from the exclude list
        self.trait_metadata = self.trait_metadata[~self.trait_metadata['rfid'].isin(exclude_rfids)]

        # drop any groups that aren't requested
        if self.n_requested_males == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata['sex'] != 'M']
        if self.n_requested_females == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata['sex'] != 'F']
        if self.n_requested_high == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata[f'{self.trait}_group'] != 'high']
        if self.n_requested_low == 0:
            self.trait_metadata = self.trait_metadata[self.trait_metadata[f'{self.trait}_group'] != 'low']

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


    def assign_rattaca(self, rfids_to_assign, rattaca_requests = None, breeders_request = None):
        '''
        Assign rats to a RATTACA project.
        
        Args:
            rfids_to_assign: List of RFIDs to assign, ordered [min_rat, max_rat].
            
            rattaca_requests: A list of all other open RATTACA request objects 
                (not including this specific request).
            
            breeders_request: HSWBreeders request, to which remaining animals 
                will be assigned by priority if needed
        '''

        if not isinstance(rfids_to_assign, list) or len(rfids_to_assign) != 2:
            raise TypeError('rfids_to_assign must be a list with two RFID elements')
        
        all_rattaca_requests = [self] + rattaca_requests if rattaca_requests is not None else self

        # validate that both rats are available
        for rfid in rfids_to_assign:
            if rfid not in self.available_rfids:
                print(f'WARNING: RFID {rfid} is not in available_rfids, skipping assignment')
                return
            if rfid not in self.trait_metadata['rfid'].values:
                print(f'ERROR: RFID {rfid} not found in trait_metadata, skipping assignment')
                return

            # make sure rats are available
        available_rfids = set(self.available_rfids)
        for rfid in rfids_to_assign:
            if rfid not in available_rfids:
                print(f'RFID {rfid} is not available for assignment')
                
        # get rat metadata
        min_rat = rfids_to_assign[0]
        max_rat = rfids_to_assign[1]
        min_rat_data = self._rfid_metadata(min_rat)
        max_rat_data = self._rfid_metadata(max_rat)
        
        min_rat_sex = min_rat_data['sex']
        min_rat_fam = min_rat_data['fam']
        min_rat_group = min_rat_data['group']
        min_rat_pred = min_rat_data['pred']
        
        max_rat_sex = max_rat_data['sex']
        max_rat_fam = max_rat_data['fam']
        max_rat_group = max_rat_data['group']
        max_rat_pred = max_rat_data['pred']
        
        n_min_rat_sibs = len(min_rat_data['sibs'])
        n_max_rat_sibs = len(max_rat_data['sibs'])

        # count the number of remaining siblings that could be assigned before
        # assigning the current rats
        remaining_min_rat_sibs_allowed = \
            self._get_n_remaining_per_fam(sex = min_rat_sex, fam = min_rat_fam)
        remaining_max_rat_sibs_allowed = \
            self._get_n_remaining_per_fam(sex = max_rat_sex, fam = max_rat_fam)

        if remaining_min_rat_sibs_allowed == 0:
            message = (
                f'RFID {rfid} cannot be assigned. \n'
                f'Breederpair {min_rat_fam} has already contributed {self.max_per_sex} {min_rat_sex} animals to {self.request_name} \n')
            print(message)
            exit
        if remaining_max_rat_sibs_allowed == 0:
            message = (
                f'RFID {rfid} cannot be assigned. \n'
                f'Breederpair {max_rat_fam} has already contributed {self.max_per_sex} {max_rat_sex} animals to {self.request_name} \n')
            print(message)
            exit
        
        # assign the high rat
        if max_rat_sex == 'M':
            if not self.is_satisfied_rattaca('M', 'high'):
                self.assigned_males_high[max_rat] = (max_rat_sex, max_rat_pred, max_rat_group)
                self.assigned_fams['M_high'][max_rat_fam] = (max_rat, max_rat_sex)
                self.remove([max_rat])
            else:
                print(f'Cannot assign {max_rat}: Male/high group already satisfied')
        elif max_rat_sex == 'F':
            if not self.is_satisfied_rattaca('F', 'high'):
                self.assigned_females_high[max_rat] = (max_rat_sex, max_rat_pred, max_rat_group)
                self.assigned_fams['F_high'][max_rat_fam] = (max_rat, max_rat_sex)
                self.remove([max_rat])
            else:
                print(f'Cannot assign {max_rat}: Female/high group already satisfied')
        
        # check if remaining siblings should be prioritized for HSW breeders
        if breeders_request is not None:
            if n_max_rat_sibs == 1:
                breeder_sib = max_rat_data['sibs']
                if max_rat_fam not in breeders_request.assigned_fams[max_rat_sex]:
                    breeders_request.assign_hsw_breeders(
                        rfids_to_assign = breeder_sib,
                        non_breeder_requests = all_rattaca_requests)

        # remove siblings from availability
        self._update_available_rats(by='fam')

        # assign the low rat
        if min_rat_sex == 'M':
            if not self.is_satisfied_rattaca('M', 'low'):
                self.assigned_males_low[min_rat] = (min_rat_sex, min_rat_pred, min_rat_group)
                self.assigned_fams['M_low'][min_rat_fam] = (min_rat, min_rat_sex)
                self.remove([min_rat])
            else:
                print(f'Cannot assign {min_rat}: Male/low group already satisfied')
        elif min_rat_sex == 'F':
            if not self.is_satisfied_rattaca('F', 'low'):
                self.assigned_females_low[min_rat] = (min_rat_sex, min_rat_pred, min_rat_group)
                self.assigned_fams['F_low'][min_rat_fam] = (min_rat, min_rat_sex)
                self.remove([min_rat])
            else:
                print(f'Cannot assign {min_rat}: Female/low group already satisfied')
        
        # check if remaining siblings should be prioritized for HSW breeders
        if breeders_request is not None:
            if n_min_rat_sibs == 1:
                breeder_sib = min_rat_data['sibs']
                if min_rat_fam not in breeders_request.assigned_fams[min_rat_sex]:
                    breeders_request.assign_hsw_breeders(
                        rfids_to_assign = breeder_sib,
                        non_breeder_requests = all_rattaca_requests)

        # remove siblings from availability
        self._update_available_rats(by='fam')

        # update delta
        self.delta += max_rat_pred - min_rat_pred
        
        # update available rats list
        self._update_available_rats(by='group')
        

    def _update_available_rats(self, by=['group','fam']):
        '''Update available rats list based on satisfied groups.'''
        
        if by == 'group':
            if self.is_satisfied_rattaca('M', 'high'):
                high_male_rfids = self.trait_metadata[
                    (self.trait_metadata['sex'] == 'M') & 
                    (self.trait_metadata[f'{self.trait}_group'] == 'high')]['rfid'].tolist()
                for rfid in high_male_rfids:
                    if rfid in self.available_rfids:
                        self.available_rfids.remove(rfid)
                        
            if self.is_satisfied_rattaca('M', 'low'):
                low_male_rfids = self.trait_metadata[
                    (self.trait_metadata['sex'] == 'M') & 
                    (self.trait_metadata[f'{self.trait}_group'] == 'low')]['rfid'].tolist()
                for rfid in low_male_rfids:
                    if rfid in self.available_rfids:
                        self.available_rfids.remove(rfid)
                        
            if self.is_satisfied_rattaca('F', 'high'):
                high_female_rfids = self.trait_metadata[
                    (self.trait_metadata['sex'] == 'F') & 
                    (self.trait_metadata[f'{self.trait}_group'] == 'high')]['rfid'].tolist()
                for rfid in high_female_rfids:
                    if rfid in self.available_rfids:
                        self.available_rfids.remove(rfid)
                        
            if self.is_satisfied_rattaca('F', 'low'):
                low_female_rfids = self.trait_metadata[
                    (self.trait_metadata['sex'] == 'F') & 
                    (self.trait_metadata[f'{self.trait}_group'] == 'low')]['rfid'].tolist()
                for rfid in low_female_rfids:
                    if rfid in self.available_rfids:
                        self.available_rfids.remove(rfid)
        
        # check availability by family x sex
        if by == 'fam':
            all_fams = self.trait_metadata['breederpair'].unique().tolist()
            for fam in all_fams:
                fam_df = self.trait_metadata[self.trait_metadata['breederpair']==fam]
                fam_males = fam_df[fam_df['sex']=='M']['rfid'].tolist()
                fam_females = fam_df[fam_df['sex']=='F']['rfid'].tolist()

                n_avail_males = self._get_n_remaining_per_fam(sex = 'M', fam = fam)
                n_avail_females = self._get_n_remaining_per_fam(sex = 'F', fam = fam)

                if n_avail_males < 1:
                    self.remove(fam_males)
                if n_avail_females < 1:
                    self.remove(fam_females)
