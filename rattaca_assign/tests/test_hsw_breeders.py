from rattaca_assignment.request import Request
from rattaca_assignment.request_hsw_breeders import HSWBreeders
from rattaca_assignment.request_rattaca import RATTACA
from rattaca_assignment.main import parse_args
import unittest
import argparse
import os
import pandas as pd
import numpy as np
import random

# set up test datasets
proj_dir = '/tscc/projects/ps-palmer/hs_rats/rattaca/assignment/rattaca_assignment.git'
req_dir = os.path.join(proj_dir, 'proj_requests')
data_dir = os.path.join(proj_dir, 'data')
metadata = os.path.join(data_dir, 'hsw_gen100_colony.csv')
preds = os.path.join(data_dir, 'rattaca_gen100_preds.csv')
exclude = os.path.join(data_dir, 'rattaca_gen100_exclude.csv')
trait_1 = os.path.join(req_dir, 'rattaca', 'trait_1.json')
trait_2 = os.path.join(req_dir, 'rattaca', 'trait_2.json')
trait_3 = os.path.join(req_dir, 'rattaca', 'trait_3.json')
hsw_breeders = os.path.join(req_dir, 'hsw_breeders', 'hsw_breeders.json')

# set up input arguments
input_args = ['-c', metadata, '-p', preds, '-r', hsw_breeders, 
    '-e', exclude]
args_excl = parse_args(input_args) # uses the list of RFIDs to exclude
args_all = parse_args(input_args[:-2]) # all rats

exclude_rfids = pd.read_csv(exclude, header = None).astype(str)[0].\
    values.tolist()

# RFIDs to test manual assignments
# (three females from breeder pair 16)
manual_rfids = ['933000320978183','933000320978184']
remaining_sib = '933000320978186'

# singletons - males or females w/o siblings of the same sex
single_males = ['933000320978204','933000320978201']
single_females = ['933000320978151','933000320946456']
singletons = single_males + single_females

# breederpairs with singleton males or females
single_m_fams = [4, 19]
single_f_fams = [46, 58]

# tests for the HSWBreeders request subclass
class TestHSWBreeders(unittest.TestCase):

    # set up the test fixture
    def setUp(self):

        # create Request instances for testing
        self.request1 = RATTACA(trait_1, args_all)
        self.hsw_request = HSWBreeders(hsw_breeders, args_all)
        self.all_requests = [self.request1, self.hsw_request]

        # read in data from arguments
        self.preds_df = pd.read_csv(preds, dtype = {'rfid': str})

        colony_df = pd.read_csv(metadata, dtype = {'rfid': str, 
            'accessid': int})
        self.exclude = pd.read_csv(exclude, header = None).astype(str)[0].values

        # identify dead rats
        dead_strs = ['dead', 'die', 'death', 'euth', 'eauth', 'kill']
        dead_search = '|'.join(dead_strs)
        dead_rats = colony_df[colony_df['comments']\
            .str.contains(dead_search, case=False, na=False)]['rfid'].tolist()
        
        # identify rats with unknown sex
        ambig_sex = colony_df[colony_df['comments']\
            .str.contains('unsure sex', case=False, na=False)]['rfid'].tolist()

        # drop dead rats, rats w/ unknown sex
        drop_rats = dead_rats + ambig_sex
        colony_df = colony_df[~colony_df['rfid'].isin(drop_rats)]

        gtyped_rfids = self.preds_df['rfid'].tolist()
        colony_df['gtyped'] = colony_df['rfid'].isin(gtyped_rfids).astype(int)

        self.colony_df = colony_df
        self.excl_df = colony_df[~colony_df['rfid'].isin(exclude_rfids)]

        self.exclude = pd.read_csv(exclude, header = None).astype(str)[0].values

        # set up metadata for each request
        trait_metadata1 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_1']], on='rfid', how='left')\
            .sort_values(by='trait_1', axis=0, ascending=False)\
            .dropna(subset='trait_1', ignore_index=True)

        excl_metadata1 = trait_metadata1[~trait_metadata1['rfid'].isin(self.exclude)]
        self.trait_md_excl = excl_metadata1

    # # test that counter properties properly track assigned vs available rats
    # def test_counters(self):

    #     req = self.hsw_request
    #     print(req.__dict__.keys())
    #     req.assign_manual(manual_rfids)

    def test_prioritize_breeders(self):
        
        # first, ensure singletons do not start off assigned
        for rfid in singletons:
            self.assertNotIn(rfid, self.hsw_request.assigned_rats['assigned_total'])

        # then ensure that all singletons are prioritized for assignment
        priority_breeders = self.hsw_request.prioritize_breeders()
        self.assertSetEqual(set(priority_breeders), set(singletons))

    def test_update_breeders(self):

        # manually assign 2 of 3 females from B16
        # first, ensure RFIDs are not assigned or prioritized to begin with
        for rfid in manual_rfids:     
            self.assertIn(rfid, self.request1.available_rfids)
            self.assertIn(rfid, self.hsw_request.available_rfids)
            self.assertNotIn(rfid, self.request1.assigned_rats['assigned_total'])
        self.assertIn(remaining_sib, self.hsw_request.available_rats)
        self.assertNotIn(remaining_sib, self.hsw_request.assigned_rats['assigned_total'])

        # then manually assign RFIDs to the RATTACA request
        self.request1.assign_manual(manual_rfids)
        self.hsw_request.remove(manual_rfids, [self.request1])

        # update the list of RFIDs now available for assignment
        new_priority = self.hsw_request.update_breeders([self.request1])
        
        # ensure the remaining sibling is now set as a priority breeder
        self.assertIn(remaining_sib, new_priority)

    def test_assign_hsw_breeders(self):
        
        # manually assign 2 of 3 females from B16
        # first, ensure RFIDs are not assigned to begin with
        for rfid in manual_rfids:     
            self.assertIn(rfid, self.request1.available_rfids)
            self.assertIn(rfid, self.hsw_request.available_rfids)
            self.assertNotIn(rfid, self.request1.assigned_rats['assigned_total'])
            self.assertNotIn(rfid, self.hsw_request.assigned_rats['assigned_total'])
        self.assertIn(remaining_sib, self.hsw_request.available_rats)
        self.assertNotIn(remaining_sib, self.hsw_request.assigned_rats['assigned_total'])

        # then manually assign RFIDs to the RATTACA request
        self.request1.assign_manual(manual_rfids)
        self.hsw_request.remove(manual_rfids, [self.request1])

        # check that assigned RFIDs are in fact assigned
        # and removed from availability
        for rfid in manual_rfids:
            self.assertIn(rfid, self.request1.assigned_rats['assigned_total'])
            self.assertNotIn(rfid, self.request1.available_rfids)
            self.assertNotIn(rfid, self.hsw_request.assigned_rats['assigned_total'])
            self.assertNotIn(rfid, self.hsw_request.available_rfids)

        # then assign breeders
        self.hsw_request.assign_hsw_breeders([self.request1])
        self.assertIn(remaining_sib, self.hsw_request.assigned_rats['assigned_total'])
        self.assertIn(remaining_sib, self.hsw_request.assigned_rats['assigned_females'])
        self.assertIn(remaining_sib, self.hsw_request.assigned_females)


    def test_assign_remainder(self):

        req = self.hsw_request

        # ensure the RFIDs are not assigned to begin with
        for rfid in singletons:
            self.assertIn(rfid, req.available_rfids)
            if rfid in single_males:
                self.assertNotIn(rfid, req.assigned_males)
            else:
                self.assertNotIn(rfid, req.assigned_females)

        # assign singletons by priority, confirm assignment
        req.assign(req.prioritize_breeders(),[])

        for rfid in singletons:
            self.assertNotIn(rfid, req.available_rfids)
            if rfid in single_males:
                self.assertIn(rfid, req.assigned_males)
            else:
                self.assertIn(rfid, req.assigned_females)

        # then assign remaining breeders
        req.assign_remainder()

        # check the number of assignments against
        # the number of breederpairs with M/F to assign
        n_male_fams = len(req.fams_with_males)
        n_female_fams = len(req.fams_with_females)

        n_males_assigned = len(req.assigned_males)
        n_females_assigned = len(req.assigned_females)

        self.assertEqual(n_male_fams, n_males_assigned)
        self.assertEqual(n_female_fams, n_females_assigned)

        # ensure the max possible number of individuals has been assigned
        n_assigned = n_males_assigned + n_females_assigned
        n_fams = len(req.all_fams)
        n_missing_m_fams = len(req.fams_without_males)
        n_missing_f_fams = len(req.fams_without_females)
        n_possible_assignments = 2*n_fams - n_missing_m_fams - n_missing_f_fams

        self.assertEqual(n_assigned, n_possible_assignments)
        
        # ensure that all families are represented where possible
        assigned_m_fams = list(req.assigned_breederpairs\
            ['assigned_fams_m'].keys())
        assigned_f_fams = list(req.assigned_breederpairs\
            ['assigned_fams_f'].keys())
        assigned_fams = set(assigned_m_fams + assigned_f_fams)
        all_fams = set(req.all_fams)

        self.assertSetEqual(assigned_fams, all_fams)

    # test the option to fill a breeder request
    def test_fill_request(self):
        
        req = self.hsw_request
        req.assign(remaining_requests = [], fill_request = True)

        # check the number of assignments against
        # the number of breederpairs with M/F to assign
        n_male_fams = len(req.fams_with_males)
        n_female_fams = len(req.fams_with_females)

        n_males_assigned = len(req.assigned_males)
        n_females_assigned = len(req.assigned_females)

        self.assertEqual(n_male_fams, n_males_assigned)
        self.assertEqual(n_female_fams, n_females_assigned)

        # ensure the max possible number of individuals has been assigned
        n_assigned = n_males_assigned + n_females_assigned
        n_fams = len(req.all_fams)
        n_missing_m_fams = len(req.fams_without_males)
        n_missing_f_fams = len(req.fams_without_females)
        n_possible_assignments = 2*n_fams - n_missing_m_fams - n_missing_f_fams

        self.assertEqual(n_assigned, n_possible_assignments)
        
        # ensure that all families are represented where possible
        assigned_m_fams = list(req.assigned_breederpairs\
            ['assigned_fams_m'].keys())
        assigned_f_fams = list(req.assigned_breederpairs\
            ['assigned_fams_f'].keys())
        assigned_fams = set(assigned_m_fams + assigned_f_fams)
        all_fams = set(req.all_fams)

        self.assertSetEqual(assigned_fams, all_fams)

    def test_is_satisfied(self):

        req = self.hsw_request

        # first, ensure the request is not satisfied to begin with
        self.assertEqual(req.is_satisfied(), False)
        self.assertEqual(req.is_satisfied(sex = 'M'), False)
        self.assertEqual(req.is_satisfied(sex = 'F'), False)
        self.assertEqual(req.is_satisfied(sex = 'M', family = 1), False)
        self.assertEqual(req.is_satisfied(sex = 'F', family = 1), False)

        # next, assign breeders by priority
        req.assign(req.prioritize_breeders(),[])

        assigned_m_fams = list(req.assigned_breederpairs\
            ['assigned_fams_m'].keys())
        assigned_f_fams = list(req.assigned_breederpairs\
            ['assigned_fams_f'].keys())

        # test that families with singletons are satisfied
        for rfid in single_males:
            self.assertIn(rfid, req.assigned_males)
        for rfid in single_females:
            self.assertIn(rfid, req.assigned_females)

        for fam in single_m_fams:
            self.assertIn(fam, assigned_m_fams)
        for fam in single_f_fams:
            self.assertIn(fam, assigned_f_fams)

        self.assertEqual(req.is_satisfied(sex='M', family=4), True)
        self.assertEqual(req.is_satisfied(sex='M', family=19), True)
        self.assertEqual(req.is_satisfied(sex='F', family=46), True)
        self.assertEqual(req.is_satisfied(sex='F', family=58), True)

        # assign remaining breeders
        req.assign_remainder()

        # the request cannot be fully satisfied with one M/F per breederpair
        # but check that individual family assignments are satisfied
        self.assertEqual(req.is_satisfied(sex='F', family=4), True)
        self.assertEqual(req.is_satisfied(sex='F', family=19), True)
        self.assertEqual(req.is_satisfied(sex='M', family=46), True)
        self.assertEqual(req.is_satisfied(sex='M', family=58), True)
        self.assertEqual(req.is_satisfied(sex='F', family=1), False)
        self.assertEqual(req.is_satisfied(sex='F', family=2), True)
        self.assertEqual(req.is_satisfied(sex='M', family=2), True)
    
        
    def test_count_leftovers(self):

        rfid = single_males[0]
        req = self.hsw_request

        # first, ensure the RFID is available, 
        # and no IDs are left over for assignment to begin with
        self.assertIn(rfid, req.available_rfids)
        self.assertNotIn(rfid, req.assigned_rats)
        self.assertEqual(req.count_leftovers()['M'], req.n_requested_males)
        self.assertEqual(req.count_leftovers()['F'], req.n_requested_females)

        # then fill the request
        req.assign_manual([rfid])
        self.assertNotIn(rfid, req.available_rfids)
        self.assertIn(rfid, req.assigned_males)
        self.assertEqual(req.count_leftovers()['M'], req.n_requested_males - 1)
        req.assign(req.prioritize_breeders(), remaining_requests=[])
        self.assertEqual(req.count_leftovers()['M'], req.n_requested_males - 2)
        self.assertEqual(req.count_leftovers()['F'], req.n_requested_females - 2)
        req.assign(remaining_requests = [], fill_request = True)
        self.assertEqual(req.count_leftovers()['M'], 6)
        self.assertEqual(req.count_leftovers()['F'], 8)
        self.assertEqual(req.count_leftovers()['total'], 14)


    def test_assign_manual(self):

        req = self.hsw_request
        sib_1 = manual_rfids[0]
        sib_2 = remaining_sib

        # first, ensure no RFIDs are assigned
        self.assertIn(sib_1, req.available_rfids)
        self.assertIn(sib_1, req.available_females)
        self.assertNotIn(sib_1, req.assigned_females)
        self.assertIn(sib_2, req.available_rfids)
        self.assertIn(sib_2, req.available_females)
        self.assertNotIn(sib_2, req.assigned_females)
        self.assertIn(16, req.available_breederpairs['available_fams_f'])
        self.assertEqual(req.is_satisfied(sex='F', family=16), False)

        # next, assign one F from breederpair 16
        req.assign_manual_hsw_breeders([sib_1])
        self.assertIn(sib_1, req.assigned_females)
        self.assertNotIn(sib_1, req.available_rfids)
        self.assertNotIn(sib_1, req.available_females)
        self.assertNotIn(16, req.available_breederpairs['available_fams_f'])
        self.assertEqual(req.is_satisfied(sex='F', family=16), True)

        # try assigning a sibling
        self.assertNotIn(sib_2, req.available_females)
        self.assertNotIn(sib_2, req.assigned_females)
        self.assertRaises(ValueError,
            req.assign_manual_hsw_breeders,[sib_2], override=False)
        req.assign_manual([sib_2], override=True)
        self.assertNotIn(sib_2, req.available_females)
        self.assertIn(sib_2, req.assigned_females)
        self.assertNotIn(remaining_sib, req.available_females)


    def test_remove(self):

        req = self.hsw_request
        sib_1 = manual_rfids[0]
        sib_2 = remaining_sib

        self.assertIn(sib_1, req.available_rfids)
        self.assertIn(sib_1, req.available_females)
        self.assertNotIn(sib_1, req.assigned_females)
        self.assertIn(16, req.available_breederpairs['available_fams_f'])
        self.assertEqual(req.is_satisfied(sex='F', family=16), False)

        # manually remove the RFID,
        # check it has been removed without assigning or excluding sibs from assignment
        req.remove([sib_1], [])
        self.assertNotIn(sib_1, req.available_rfids)
        self.assertNotIn(sib_1, req.available_females)
        self.assertNotIn(sib_1, req.assigned_females)
        self.assertIn(16, req.available_breederpairs['available_fams_f'])
        self.assertEqual(req.is_satisfied(sex='F', family=16), False)


if __name__ == '__main__':
    unittest.main()