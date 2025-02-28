### dependencies ###
import json
import random
from rattaca_assignment.request import Request

### Request subclass for HS West breeders
class HSWBreeders(Request):

    def __init__(self, request_file, args):
                
        # inherit attributes from the parent Request class
        super().__init__(request_file, args)

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

        # initialize lists of breeder pairs with rats available for assignment
        # as breeders
        self.available_fams_m = {}
        self.available_fams_f = {}
        
        self.all_fams = self.colony_df['breederpair'].unique()
        fams_with_males = self.colony_df.groupby('breederpair').\
            filter(lambda df: 'M' in df['sex'].values)
        fams_with_females = self.colony_df.groupby('breederpair').\
            filter(lambda df: 'F' in df['sex'].values)

        self.fams_with_males = fams_with_males['breederpair'].unique()
        self.fams_with_females = fams_with_females['breederpair'].unique()

        self.fams_without_males = \
            list(set(self.all_fams) - set(self.fams_with_males))
        self.fams_without_females = \
            list(set(self.all_fams) - set(self.fams_with_females))


    def update_breeders(self, non_breeder_requests):

        if isinstance(non_breeder_requests, list):
            pass
        else:
            raise TypeError('non_breeder_requests must be a list of Request classes')

        use_requests = [req for req in non_breeder_requests if not req.is_satisfied()]
        df = self.trait_metadata

        # list all rats that have been assigned to other projects
        assigned_males = []
        assigned_females = []
        for request in use_requests:
            request_males = \
                list(request.assigned_rats['assigned_males'].keys())
            request_females = \
                list(request.assigned_rats['assigned_females'].keys())
            assigned_males.extend(request_males)
            assigned_females.extend(request_females)

        male_breeders = list(self.assigned_males.keys())
        female_breeders = list(self.assigned_females.keys())
        
        filled_breeder_fams_m = self.assigned_breederpairs['assigned_fams_m']
        filled_fam_rfids_m = [v[1] for v in filled_breeder_fams_m.values()] 
        filled_breeder_fams_f = self.assigned_breederpairs['assigned_fams_f']
        filled_fam_rfids_f = [v[1] for v in filled_breeder_fams_f.values()] 

        # list all rats not available for assignment: 
        # breeders + breeder siblings + rats assigned to other projects
        unavail_males = \
            assigned_males + male_breeders + filled_fam_rfids_m
        unavail_females = \
            assigned_females + female_breeders + filled_fam_rfids_f
        unavail_rfids = unavail_males + unavail_females

        # update the list of rats available for assignment to breeders
        self.available_rfids = \
            list(set(self.available_rfids) - set(unavail_rfids))

        # check if any rats are now priority for assignment to breeders
        new_priority_breeders = self.prioritize_breeders()

        return new_priority_breeders


    # function to assign breeders a priori
    def prioritize_breeders(self):
        '''Assigns singleton offspring to HSW breeders by priority 
        prior to project assignments.
        
        Identifies all rats that are the only M or F from their litter,
        assigns them automatically to HSW breeders to ensure all breeder
        pairs have one M and one F represented in the next HSW generation.

        Args:
            df: The 'trait_metadata' pandas dataframe for HSW breeders
        
        Returns:
            A list of RFIDs that must be assigned to HSW breeders.
        '''

        # only examine rats that are available for assignment
        df = self.trait_metadata
        df = df[df['rfid'].isin(self.available_rfids)]
        m_breeders = []
        f_breeders = []

        for pair in df['breederpair'].unique():

            litter_males = df[(df['breederpair'] == pair) \
                & (df['sex'] == 'M')]['rfid'].tolist()
            litter_females = df[(df['breederpair'] == pair) \
                & (df['sex'] == 'F')]['rfid'].tolist()
            
            # any rat that is the single M or F from its litter 
            # is automatically assigned to HS West breeders
            if len(litter_males) == 1:
                assign_male = litter_males
                m_breeders.extend(assign_male)
            
            if len(litter_females) == 1:
                assign_female = litter_females
                f_breeders.extend(assign_female)
        
        breeders = m_breeders + f_breeders
        return(breeders)


    def assign_hsw_breeders(self, non_breeder_requests=None, fill_request=False):

        # assign breeders by priority
        breeders_to_assign = self.prioritize_breeders()
        breeders_to_assign = {k: v for k, v in self.available_rats.items() \
            if k in breeders_to_assign}
        males_to_assign = {k: v for k, v in breeders_to_assign.items() \
            if v[0] == 'M'}
        females_to_assign = {k: v for k, v in breeders_to_assign.items() \
            if v[0] == 'F'}

        # if the request is not yet satisfied
        # assign rats and remove them from availability
        successfully_assigned = []
        for assigned_breeder in breeders_to_assign.items():
            breeder_rfid = assigned_breeder[0]
            breeder_sex = assigned_breeder[1][0]
            breeder_fam = assigned_breeder[1][1]
            
            # assign the rat and remove it from availability
            if breeder_sex == 'M':
                assign_to = self.assigned_males
            elif breeder_sex == 'F':
                assign_to = self.assigned_females

            if self.is_satisfied_hsw_breeders(sex = breeder_sex, family = breeder_fam):
                message = (f'RFID {breeder_rfid} could not be assigned:'
                f'{breeder_sex} assignments for breeder pair {breeder_fam}'
                'already satisfied')
                print(message)
            else:
                assign_to[breeder_rfid] = breeders_to_assign[breeder_rfid]
                successfully_assigned.append(breeder_rfid)
                # print(f'Assigned {breeder_sex} breeder {breeder_rfid} from pair {breeder_fam} by priority')
        
        self.remove(
            rats_to_remove = successfully_assigned, 
            remaining_requests = non_breeder_requests)

        # fill all remaining breeder assignments if desired
        if fill_request is True:
            self.assign_remainder()
    

    # fill all remaining breeder assignments
    def assign_remainder(self):
        # print('ASSIGN_REMAINDER')
        final_males = {}
        final_females = {}
        avail_males = self.available_breederpairs['available_fams_m']
        avail_females = self.available_breederpairs['available_fams_f']
        
        # randomly sample one male from each available breeder pair
        for fam, males in avail_males.items():
            if males:
                final_males[fam] = \
                    random.choice([male[0] for male in males])
        for fam, females in avail_females.items():
            if females:
                final_females[fam] = \
                    random.choice([female[0] for female in females])

        final_males = list(final_males.values())
        final_females = list(final_females.values())
        all_breeders = final_males + final_females
        final_breeders = {rfid: self.available_rats[rfid] \
            for rfid in all_breeders if rfid in self.available_rfids}
        
        # if the request is not yet satisfied
        # assign rats and remove them from availability
        successfully_assigned = []
        for assigned_breeder in final_breeders.items():
            breeder_rfid = assigned_breeder[0]
            breeder_sex = assigned_breeder[1][0]
            breeder_fam = assigned_breeder[1][1]
            
            # assign the rat and remove it from availability
            if breeder_sex == 'M':
                assign_to = self.assigned_males
            elif breeder_sex == 'F':
                assign_to = self.assigned_females

            if self.is_satisfied(sex=breeder_sex, family=breeder_fam):
                message = (
                    f'RFID {breeder_rfid} could not be assigned: {breeder_sex} \
                        assignments for breeder pair {breeder_fam} already satisfied'
                )
                print(message)

            else:
                assign_to[breeder_rfid] = final_breeders[breeder_rfid]
                successfully_assigned.append(breeder_rfid)
                # print(f'Assigned {breeder_sex} breeder: {breeder_rfid} from pair {breeder_fam} to fill remainder')
            
        self.remove(
            rats_to_remove = successfully_assigned, 
            remaining_requests = [])


    # function to manually assign rats to a projects as needed
    def assign_manual_hsw_breeders(self, rfids_to_assign, override = False):

        if isinstance(rfids_to_assign, list):
            pass
        else:
            raise TypeError('rfids_to_assign must be a list')
        
        # sample trait metadata for rats to assign
        to_assign_metadata = self.trait_metadata\
            [self.trait_metadata['rfid'].isin(rfids_to_assign)]
        # convert the metadata df to a dictionary 
        to_assign_dict = to_assign_metadata.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        rats_to_assign = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in to_assign_dict.items()}

        for rfid in rats_to_assign:

            rat_sex =  rats_to_assign[rfid][0]
            rat_fam =  rats_to_assign[rfid][1]

            if rat_sex == 'M':
                assign_to = self.assigned_males
                already_assigned_sibs = \
                [rfid for rfid, (sex, fam) in self.assigned_males.items() \
                    if sex == 'M' and fam == rat_fam]
            elif rat_sex == 'F':
                assign_to = self.assigned_females
                already_assigned_sibs = \
                [rfid for rfid, (sex, fam) in self.assigned_females.items() \
                    if sex == 'F' and fam == rat_fam]

            # assign if no siblings have yet been assigned
            # for the current breederpair and sex
            if len(already_assigned_sibs) == 0:
                assign_to[rfid] = rats_to_assign[rfid]
                self.remove([rfid], remaining_requests = [])
                print(f'Manually assigned {rfid} to HS West breeders')


            # if the current RFID is from a breederpair that has already
            # assigned a sibling of the same sex to be a breeder, 
            # either prevent assignment or raise a warning            
            else:
                n_sibs = len(already_assigned_sibs)

                if override is False:
                    message = (
                        f'RFID {rfid} is not available for assignment. \n'
                        f'{n_sibs} {rat_sex} rat has already been assigned from breeder pair {rat_fam} (RFID: {already_assigned_sibs})\n'
                        f'Use override = True to force assignment of this RFID \n')
                    raise ValueError(message)

                else:
                    assign_to[rfid] = rats_to_assign[rfid]
                    self.remove([rfid], remaining_requests = [])
                    message = (
                        f'Manually assigned {rfid} to HS West breeders \n'
                        f'!!! \n'
                        f'!!! NOTE: Breeder pair {rat_fam} already had {rat_sex} sibling(s) assigned: {already_assigned_sibs})\n'
                        f'!!! Breederpair {rat_fam} now has {n_sibs + 1} {rat_sex} breeders assigned \n'
                        f'!!! \n')
                    print(message)


    def is_satisfied_hsw_breeders(self, sex=None, family=None):

        all_fams = self.colony_df['breederpair'].unique()
        fams_with_males = self.colony_df.groupby('breederpair').\
            filter(lambda df: 'M' in df['sex'].values)
        fams_with_males = fams_with_males['breederpair'].unique()
        fams_with_females = self.colony_df.groupby('breederpair').\
            filter(lambda df: 'F' in df['sex'].values)
        fams_with_females = fams_with_females['breederpair'].unique()

        if family == None:
 
            if sex is None:
                n_remaining = self.n_remaining['n_remaining_total']
            elif sex == 'M':
                n_remaining = self.n_remaining['n_remaining_males']
            elif sex == 'F':
                n_remaining = self.n_remaining['n_remaining_females']
            else:
                raise TypeError('sex must be M, F, or None')

        else:
            if family not in all_fams:
                raise ValueError(
                    f'''Could not find breeder pair {family} in the colony
                    dataframe''')

            available_m_fams = self.available_breederpairs['available_fams_m']
            available_f_fams = self.available_breederpairs['available_fams_f']

            if sex is None:
                if family in self.fams_with_males and \
                    family in self.fams_without_females:
                    if family not in available_m_fams:
                        n_remaining = 0
                    else:
                        fam_males = available_m_fams[family]
                        n_remaining = len(fam_males)
                elif family in self.fams_with_females and \
                    family in self.fams_without_males:
                    if family not in available_f_fams:
                        n_remaining = 0
                    else:
                        fam_females = available_f_fams[family]
                        n_remaining = len(fam_females)
            
            elif sex == 'M':
                if family in self.fams_without_males:
                    n_remaining = 9999
                    print(f'Breeder pair {family} has no males to assign')
                else:
                    if family not in available_m_fams:
                        n_remaining = 0
                    else: 
                        fam_males = available_m_fams[family]
                        n_remaining = len(fam_males)
            
            elif sex == 'F':
                if family in self.fams_without_females:
                    n_remaining = 9999
                    print(f'Breeder pair {family} has no females to assign')
                else:
                    if family not in available_f_fams:
                        n_remaining = 0
                    else:
                        fam_females = available_f_fams[family]
                        n_remaining = len(fam_females)

        # manually set breeders as 'satisfied' if all possible assignments 
        # have been made algorithmically
        if family is None and sex is None:
            # leftovers = self.count_leftovers()
            # if leftovers is not None:
            #     n_remaining = leftovers['total']
            self.count_leftovers()

        return(n_remaining == 0)

    def count_leftovers(self):

            n_remaining = self.n_remaining['n_remaining_total']
            n_assigned_males = len(self.assigned_males)
            n_assigned_females = len(self.assigned_females)
            n_assigned_total = n_assigned_males + n_assigned_females
            n_fams_with_males = len(self.fams_with_males)
            n_fams_with_females = len(self.fams_with_females)
            n_requested_males = self.n_requested_males
            n_requested_females = self.n_requested_females

            out = {'M': self.n_remaining['n_remaining_males'],
                   'F': self.n_remaining['n_remaining_females'],
                   'total': n_remaining}

            # if males and females have been assigned from all available breeder pairs
            if n_remaining > 0 and \
                n_assigned_males == n_fams_with_males and \
                n_assigned_females == n_fams_with_females:
                
                m_remainder = n_requested_males - n_assigned_males
                f_remainder = n_requested_females - n_assigned_females

                message = (
                    f'\n\n'
                    f'########################################################################### \n'
                    f'########################################################################### \n'
                    f'\n'
                    f'#################### HSW BREEDER ASSIGNMENT INCOMPLETE #################### \n'
                    f'\n'
                    f'All {n_assigned_total} possible breeder assignments have been made algorithmically: \n'
                    f'\t {n_assigned_males} males assigned from {n_fams_with_males} available breeder pairs \n'
                    f'\t {n_assigned_females} females assigned from {n_fams_with_females} available breeder pairs \n'
                    f'\n'
                    f'{n_requested_males} male breeders requested, {n_assigned_males} assigned. \n'
                    f'\t {m_remainder} male siblings must be assigned manually \n'
                    f'\n'
                    f'{n_requested_females} female breeders requested. {n_assigned_females} assigned. \n'
                    f'\t {f_remainder} female siblings must be assigned manually \n.'
                    f'\n'
                    f'########################################################################### \n'
                    f'########################################################################### \n\n'
                )
                print(message)
                out = {'M': m_remainder, 
                               'F': f_remainder, 
                               'total': m_remainder + f_remainder}

            return out


    # property to track rats that have been assigned to the project
    @property
    def assigned_rats(self):
        assigned_males = self.assigned_males
        assigned_females = self.assigned_females
        assigned_total = assigned_males | assigned_females
        out = {'assigned_males': assigned_males, 
               'assigned_females': assigned_females,
               'assigned_total': assigned_total}
        return out

    # property to keep track of all breeder pairs currently available
    # to sample for assignment
    @property
    def assigned_breederpairs(self):

        # update trait metadata to include only breederpairs that are currently 
        # available for assignment
        assigned_rats = self.assigned_rats['assigned_total']
        current_metadata = self.trait_metadata[
            self.trait_metadata['rfid'].isin(assigned_rats)]
        
        # convert the metadata df to a dictionary 
        # with key:breederpair and value:(rfid, sex)
        current_metadata_dict = current_metadata\
            .groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()

        # subset the dictionary to separate males and females
        assigned_fams_m = {}
        assigned_fams_f = {}

        if self.n_assigned['n_assigned_males'] == 0:
            pass
        else:
            for fam, rats in current_metadata_dict.items():
                for rfid, sex in rats:
                    if sex == 'M':
                        if fam not in assigned_fams_m:
                            assigned_fams_m[fam] = ()
                        assigned_fams_m[fam] += (rfid, sex)
        
        if self.n_assigned['n_assigned_females'] == 0:
            pass
        else:
            for fam, rats in current_metadata_dict.items():
                for rfid, sex in rats:
                    if sex == 'F':
                        if fam not in assigned_fams_f:
                            assigned_fams_f[fam] = ()
                        assigned_fams_f[fam] += (rfid, sex)
        out = {'assigned_fams_m' : assigned_fams_m, \
            'assigned_fams_f' : assigned_fams_f}
            
        return out

    # property to keep track of all breeder pairs currently available
    # to sample for assignment
    @property
    def available_breederpairs(self):

        assigned_fams_m = list(self.assigned_breederpairs\
            ['assigned_fams_m'].keys())
        assigned_fams_f = list(self.assigned_breederpairs\
            ['assigned_fams_f'].keys())

        available_fams_m = set(self.fams_with_males) - \
            set(assigned_fams_m) - set(self.fams_without_males)
        available_fams_f = set(self.fams_with_females) - \
            set(assigned_fams_f) - set(self.fams_without_females)

        current_metadata_m = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_m)]
        current_metadata_m = current_metadata_m\
            [current_metadata_m['sex'] == 'M']
        current_metadata_f = \
            self.trait_metadata[self.trait_metadata['breederpair']\
                .isin(available_fams_f)]
        current_metadata_f = current_metadata_f\
            [current_metadata_f['sex'] == 'F']

        available_fams_m = \
            current_metadata_m.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()
        available_fams_f = \
            current_metadata_f.groupby('breederpair')[['rfid', 'sex']]\
            .apply(lambda x: tuple(list(zip(x['rfid'], x['sex'])))).to_dict()

        out = {'available_fams_m': available_fams_m,
               'available_fams_f': available_fams_f}

        return out

    # property to track all non-assigned male siblings 
    # currently excluded from assignment
    @property
    def unavail_male_sibs(self):
        assigned_m_fams = \
            list(self.assigned_breederpairs['assigned_fams_m'].keys())
        m_metadata = self.trait_metadata[self.trait_metadata['sex']=='M']
        assigned_m_fams_md = \
            m_metadata[m_metadata['breederpair'].isin(assigned_m_fams)]
        current_m_sibs_md = assigned_m_fams_md[
            ~assigned_m_fams_md['rfid'].isin(self.assigned_males)]
        current_m_sibs_dict = current_m_sibs_md.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        unavail_males = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_m_sibs_dict.items()}
        
        return unavail_males

    # property to track all non-assigned female siblings 
    # currently excluded from assignment
    @property
    def unavail_female_sibs(self):
        assigned_f_fams = \
            list(self.assigned_breederpairs['assigned_fams_f'].keys())
        f_metadata = self.trait_metadata[self.trait_metadata['sex']=='F']
        assigned_f_fams_md = \
            f_metadata[f_metadata['breederpair'].isin(assigned_f_fams)]
        current_f_sibs_md = assigned_f_fams_md[
            ~assigned_f_fams_md['rfid'].isin(self.assigned_females)]
        current_f_sibs_dict = current_f_sibs_md.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        unavail_females = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_f_sibs_dict.items()}
        
        return unavail_females

    # property to track all males currently available for assignment
    @property
    def available_males(self):
        assigned_m_fams = \
            list(self.assigned_breederpairs['assigned_fams_m'].keys())
        m_metadata = self.trait_metadata[self.trait_metadata['sex']=='M']
        current_m_metadata = \
            m_metadata[~m_metadata['breederpair'].isin(assigned_m_fams)]
        current_m_md_dict = current_m_metadata.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        available_males = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_m_md_dict.items()}
        available_males = {k: v for k, v in available_males.items() \
            if k in self.available_rfids}

        return available_males
    
    # property to track all females currently available for assignment
    @property
    def available_females(self):
        assigned_f_fams = \
            list(self.assigned_breederpairs['assigned_fams_f'].keys())
        f_metadata = self.trait_metadata[self.trait_metadata['sex']=='F']
        current_f_metadata = \
            f_metadata[~f_metadata['breederpair'].isin(assigned_f_fams)]
        current_f_md_dict = current_f_metadata.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        available_females = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_f_md_dict.items()}
        available_females = {k: v for k, v in available_females.items() \
            if k in self.available_rfids}
        
        return available_females

    # property to track all rats currently available for assignment
    @property
    def available_rats(self):
        return self.available_males | self.available_females

   # counters to track the number of assigned rats
    @property
    def n_assigned(self):
        n_assigned_males = len(self.assigned_rats['assigned_males'])
        n_assigned_females = len(self.assigned_rats['assigned_females'])
        n_assigned_total = len(self.assigned_rats['assigned_total'])

        out = {'n_assigned_males': n_assigned_males,
               'n_assigned_females': n_assigned_females,
               'n_assigned_total': n_assigned_total}
        
        return out

    # counters to track the number of assignments remaining
    @property
    def n_remaining(self):

        n_remaining_males = self.n_requested_males - \
            self.n_assigned['n_assigned_males']
        n_remaining_females = self.n_requested_females - \
            self.n_assigned['n_assigned_females']
        n_remaining_total = self.n_requested_total - \
            self.n_assigned['n_assigned_total']

        out = {'n_remaining_males': n_remaining_males,
               'n_remaining_females': n_remaining_females,
               'n_remaining_total': n_remaining_total}
        
        return out
