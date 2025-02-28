### dependencies ###
import json
import random
from rattaca_assignment.request import Request

### Request subclass for randomly-assigned projects
class RandProj(Request):
    
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
            self.max_males_per_fam = req_metadat['max_males_per_fam']
            self.max_females_per_fam = req_metadat['max_females_per_fam']

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


    # function to check if a randomly-assigned project has 
    # successfully assigned all requested rats
    def is_satisfied_random(self, sex=None):

        if sex is None:
            n_remaining = self.n_remaining['n_remaining_total']
        if sex == 'M':
            n_remaining = self.n_remaining['n_remaining_males']
        elif sex == 'F':
            n_remaining = self.n_remaining['n_remaining_females']

        return(n_remaining == 0)


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

    @property
    def n_remaining_per_fam(self):

        males_df = self.trait_metadata[self.trait_metadata['sex'] == 'M']
        females_df = self.trait_metadata[self.trait_metadata['sex'] == 'F']
        total_males_per_fam = males_df\
            .groupby('breederpair')['rfid'].count().to_dict()
        total_females_per_fam = females_df\
            .groupby('breederpair')['rfid'].count().to_dict()

        unassigned_males = \
            males_df[~males_df['rfid'].isin(self.assigned_males)]
        n_unassigned_fam_males = \
            unassigned_males.groupby('breederpair')['rfid'].count().to_dict()
        
        unassigned_females = \
            females_df[~females_df['rfid'].isin(self.assigned_females)]
        n_unassigned_fam_females = \
            unassigned_females.groupby('breederpair')['rfid'].count().to_dict()

        available_males = {}
        available_females = {}

        for fam, n_unassigned in n_unassigned_fam_males.items():
            n_sibs = total_males_per_fam[fam]
            already_assigned = self.assigned_males.get(fam, 0)
            max_available = \
                min(n_unassigned, self.max_males_per_fam - already_assigned)
            available_males[fam] = max_available if max_available > 0 else 0

        for fam, n_unassigned in n_unassigned_fam_females.items():
            n_sibs = total_females_per_fam[fam]
            already_assigned = self.assigned_females.get(fam, 0)
            max_available = \
                min(n_unassigned, self.max_females_per_fam - already_assigned)
            available_females[fam] = max_available if max_available > 0 else 0

        out = {'n_available_males': available_males,
               'n_available_females': available_females}
        
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

        n_males_remaining = \
            self.n_remaining_per_fam['n_available_males']
        n_females_remaining = \
            self.n_remaining_per_fam['n_available_females']

        # list all families who still have males or females available
        # for assignment to the project
        fams_without_avail_males = \
            [k for k, v in n_males_remaining.items() if v == 0]
        fams_without_avail_females = \
            [k for k, v in n_females_remaining.items() if v == 0]

        # breeder pairs with (male) rats available to assign =
        # all pairs with males - all pairs w/o males - all pairs 
        # that have already assigned the max number of males allowed
        available_fams_m = set(self.fams_with_males) - \
            set(self.fams_without_males) - set(fams_without_avail_males)
        available_fams_f = set(self.fams_with_females) - \
            set(self.fams_without_females) - set(fams_without_avail_females)

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

    @property
    def filled_breederpairs(self):

        assigned_m_fams = self.assigned_breederpairs['assigned_fams_m']
        assigned_f_fams = self.assigned_breederpairs['assigned_fams_f']

        available_m_fams = self.available_breederpairs['available_fams_m']
        available_f_fams = self.available_breederpairs['available_fams_f']

        filled_fams_m = {}
        filled_fams_f = {}

        for fam in self.fams_with_males:
            if fam in assigned_m_fams & fam not in available_m_fams:
                filled_fams_m[fam] = assigned_m_fams[fam]
        for fam in self.fams_with_females:
            if fam in assigned_f_fams & fam not in available_f_fams:
                filled_fams_f[fam] = assigned_f_fams[fam]
        
        out = {'filled_fams_m': filled_fams_m,
               'filled_fams_f': filled_fams_f}
        
        return out


    # property to track all non-filled male siblings 
    # currently excluded from assignment
    @property
    def unavail_male_sibs(self):
        filled_m_fams = \
            list(self.filled_breederpairs['filled_fams_m'].keys())
        m_metadata = self.trait_metadata[self.trait_metadata['sex']=='M']
        filled_m_fams_md = \
            m_metadata[m_metadata['breederpair'].isin(filled_m_fams)]
        current_m_sibs_md = filled_m_fams_md[
            ~filled_m_fams_md['rfid'].isin(self.assigned_males)]
        current_m_sibs_dict = current_m_sibs_md.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        unavail_males = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_m_sibs_dict.items()}
        
        return unavail_males

    # property to track all non-filled female siblings 
    # currently excluded from assignment
    @property
    def unavail_female_sibs(self):
        filled_f_fams = \
            list(self.filled_breederpairs['filled_fams_f'].keys())
        f_metadata = self.trait_metadata[self.trait_metadata['sex']=='F']
        filled_f_fams_md = \
            f_metadata[f_metadata['breederpair'].isin(filled_f_fams)]
        current_f_sibs_md = filled_f_fams_md[
            ~filled_f_fams_md['rfid'].isin(self.assigned_females)]
        current_f_sibs_dict = current_f_sibs_md.\
            set_index('rfid')[['sex','breederpair']].to_dict(orient='index')
        unavail_females = \
            {k: (v['sex'], v['breederpair']) \
                for k, v in current_f_sibs_dict.items()}
        
        return unavail_females

    # property to track all males currently available for assignment
    @property
    def available_males(self):

        available_males = {}
        avail_m_fams = self.available_breederpairs['available_fams_m']

        for fam, males in avail_m_fams.items():
            for rfid, sex in males:
                available_males[rfid] = (sex, fam)

        return available_males
    
    # property to track all females currently available for assignment
    @property
    def available_females(self):

        available_females = {}
        avail_f_fams = self.available_breederpairs['available_fams_f']

        for fam, females in avail_f_fams.items():
            for rfid, sex in females:
                available_females[rfid] = (sex, fam)

        return available_females
        
    # property to track all rats currently available for assignment
    @property
    def available_rats(self):
        return self.available_males | self.available_females


    # update the list of available RFIDs following removal
    def update_random(self, remaining_requests):

        if isinstance(non_breeder_requests, list):
            pass
        else:
            raise TypeError('remaining_requests must be a list of Request classes')

        use_requests = [req for req in remaining_requests if not req.is_satisfied()]
        df = self.trait_metadata

        # list all rats that have been assigned to other projects
        males_assigned_elsewhere = []
        females_assigned_elsewhere = []
        for request in use_requests:
            request_males = \
                list(request.assigned_rats['assigned_males'].keys())
            request_females = \
                list(request.assigned_rats['assigned_females'].keys())
            males_assigned_elsewhere.extend(request_males)
            females_assigned_elsewhere.extend(request_females)

        assigned_males = list(self.assigned_males.keys())
        assigned_females = list(self.assigned_females.keys())
        
        filled_breeder_fams_m = self.filled_breederpairs['filled_fams_m']
        filled_fam_rfids_m = [v[1] for v in filled_breeder_fams_m.values()] 
        filled_breeder_fams_f = self.filleed_breederpairs['assigned_fams_f']
        filled_fam_rfids_f = [v[1] for v in filled_breeder_fams_f.values()] 

        # list all rats not available for assignment: 
        # breeders + breeder siblings + rats assigned to other projects
        unavail_males = \
            assigned_males + males_assigned_elsewhere + filled_fam_rfids_m
        unavail_females = \
            assigned_females + females_assigned_elsewhere + filled_fam_rfids_f
        unavail_rfids = unavail_males + unavail_females

        # update the list of rats available for assignment to breeders
        self.available_rfids = \
            list(set(self.available_rfids) - set(unavail_rfids))

    def assign_randproj(self, remaining_requests=None, fill_request=False):

        available_m_fams = self.available_breederpairs['available_fams_m']
        available_f_fams = self.available_breederpairs['available_fams_f']
        
        # randomly sample one available breeder pair per sex
        random_m_fam = random.choice([fam for fam in available_m_fams.items()])
        random_f_fam = random.choice([fam for fam in available_f_fams.items()])

        sampled_m_fam = random_m_fam[0]
        fam_males = random_m_fam[1]
        sampled_f_fam = random_f_fam[0]
        fam_females = random_f_fam[1]

        # randomly sample one rat from each family
        male_to_assign = \
            random.choice([male for male in fam_males])
        female_to_assign = \
            random.choice([female for female in fam_females])
        rats_to_assign = [male_to_assign, female_to_assign]
        print('rats_to_assign:', rats_to_assign)

        # if the request is not yet satisfied
        # assign rats and remove them from availability
        successfully_assigned = []
        for assigned_rat in rats_to_assign:
            rat_rfid = assigned_rat[0]
            rat_sex = assigned_rat[1]
            
            # assign the rat and remove it from availability
            if breeder_sex == 'M':
                assign_to = self.assigned_males
                rat_fam = sampled_m_fam
            elif breeder_sex == 'F':
                assign_to = self.assigned_females
                rat_fam = sampled_f_fam

            if self.is_satisfied_hsw_breeders(sex = rat_sex, family = rat_fam):
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
            remaining_requests = remaining_requests)

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
        all_rats = final_males + final_females
        final_rats = {rfid: self.available_rats[rfid] \
            for rfid in all_rats if rfid in self.available_rfids}
        
        # if the request is not yet satisfied
        # assign rats and remove them from availability
        successfully_assigned = []
        for assigned_rat in final_rats.items():
            rat_rfid = assigned_rat[0]
            rat_sex = assigned_rat[1][0]
            rat_fam = assigned_rat[1][1]
            
            # assign the rat and remove it from availability
            if rat_sex == 'M':
                assign_to = self.assigned_males
            elif rat_sex == 'F':
                assign_to = self.assigned_females

            if self.is_satisfied(sex=rat_sex, family=rat_fam):
                message = (
                    f'RFID {rat_rfid} could not be assigned: {rat_sex} \
                        assignments for breeder pair {rat_fam} already satisfied'
                )
                print(message)

            else:
                assign_to[rat_rfid] = final_rats[rat_rfid]
                successfully_assigned.append(rat_rfid)
                # print(f'Assigned {breeder_sex} breeder: {breeder_rfid} from pair {breeder_fam} to fill remainder')
            
        self.remove(
            rats_to_remove = successfully_assigned, 
            remaining_requests = [])

