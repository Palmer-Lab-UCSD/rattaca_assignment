from rattaca_assignment.request import Request
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

# set up input arguments
input_args = ['-c', metadata, '-p', preds, '-r', trait_1, trait_2, 
              trait_3, '-e', exclude]
args_excl = parse_args(input_args) # uses the list of RFIDs to exclude
args_all = parse_args(input_args[:-2]) # all rats

exclude_rfids = pd.read_csv(exclude, header = None).astype(str)[0].\
    values.tolist()

### FINISH: test that group assignment works
def get_trait_groups(df, trait, n_groups=2):

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

# tests for argument parsing
class TestParser(unittest.TestCase):
    
    # test that args elements are read in correctly
    def test_parse_args_simple(self):
        args = parse_args(['-c', 'colony_df_file',
                              '-p', 'preds_file',
                              '-r','first_request_file','second_request_file'])
        self.assertEqual(args.colony_dataframe[0], 'colony_df_file')
        self.assertEqual(args.predictions[0], 'preds_file')
        self.assertEqual(args.requests[0], 'first_request_file')
        self.assertEqual(args.requests[1], 'second_request_file')    
    

# tests for the Request parent class
class TestRequest(unittest.TestCase):

    # set up the test fixture
    def setUp(self):

        # create Request instances for testing
        self.request1 = Request(trait_1, args_all) # all rats
        self.request2 = Request(trait_2, args_all)
        self.request3 = Request(trait_3, args_all)
        self.excl_req1 = Request(trait_1, args_excl) # uses the exclude file
        self.excl_req2 = Request(trait_2, args_excl) # uses the exclude file
        self.excl_req3 = Request(trait_3, args_excl) # uses the exclude file

        self.full_req_list = [self.request1, self.request2, self.request3]
        self.excl_req_list = [self.excl_req1, self.excl_req2, self.excl_req3]
        self.all_requests = self.full_req_list + self.excl_req_list

        # read in colony data from arguments
        preds_df = pd.read_csv(preds, dtype = {'rfid': str})
        gtyped_rfids = preds_df['rfid'].tolist()
        colony_df = pd.read_csv(metadata, dtype = {'rfid': str, 
            'accessid': int})
        self.exclude = pd.read_csv(exclude, header = None).astype(str)[0].values

        # identify rats that have been genotyped
        colony_df['gtyped'] = colony_df['rfid'].isin(gtyped_rfids).astype(int)
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

        # drop from exclude list
        self.colony_df = colony_df[~colony_df['rfid'].isin(drop_rats)]
        self.excl_df = self.colony_df[~self.colony_df['rfid'].isin(exclude_rfids)]

    # test that the colony df is read in properly, with excluded RFIDs dropped
    def test_colony_df(self):
        
        full_df = self.colony_df
        excl_df = full_df[~full_df['rfid'].isin(exclude_rfids)]
        
        for req in self.full_req_list:
            pd.testing.assert_frame_equal(req.colony_df, full_df)

        for req in self.excl_req_list:
            pd.testing.assert_frame_equal(req.colony_df, full_df)        

    # def test_assign(self):
    #     with self.assertRaises(NotImplementedError):
    #         self.request1.update()

#     def test_is_satistfied(self):
#         with self.assertRaises(NotImplementedError):
#             self.request1.is_satisfied()

    # # tidy up the test fixture after testing
    # def tearDown(self):
    #     self.request1.dispose()



# tests for the RATTACA request subclass
class TestRATTACA(unittest.TestCase):

    # set up the test fixture
    def setUp(self):

        # create Request instances for testing
        self.request1 = RATTACA(trait_1, args_all) # all rats
        self.request2 = RATTACA(trait_2, args_all)
        self.request3 = RATTACA(trait_3, args_all)
        self.excl_req1 = RATTACA(trait_1, args_excl) # uses the exclude file
        self.excl_req2 = RATTACA(trait_2, args_excl) # uses the exclude file
        self.excl_req3 = RATTACA(trait_3, args_excl) # uses the exclude file

        self.full_req_list = [self.request1, self.request2, self.request3]
        self.excl_req_list = [self.excl_req1, self.excl_req2, self.excl_req3]
        self.all_requests = self.full_req_list + self.excl_req_list

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
        trait_metadata2 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_2']], on='rfid', how='left')\
            .sort_values(by='trait_2', axis=0, ascending=False)\
            .dropna(subset='trait_2', ignore_index=True)
        trait_metadata3 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_3']], on='rfid', how='left')\
            .sort_values(by='trait_3', axis=0, ascending=False)\
            .dropna(subset='trait_3', ignore_index=True)
        
        # add group assignments to each request's metadata
        trait_metadata1['trait_1_2group'] = \
            get_trait_groups(trait_metadata1, 'trait_1', n_groups=2)
        trait_metadata2['trait_2_2group'] = \
            get_trait_groups(trait_metadata2, 'trait_2', n_groups=2)
        # drop unassigned groups
        trait_metadata2 = trait_metadata2[trait_metadata2['sex'] != 'M']
        trait_metadata3['trait_3_2group'] = \
            get_trait_groups(trait_metadata3, 'trait_3', n_groups=2)
        self.trait_md = [trait_metadata1, trait_metadata2, trait_metadata3]

        excl_metadata1 = trait_metadata1[~trait_metadata1['rfid'].isin(self.exclude)]
        excl_metadata2 = trait_metadata2[~trait_metadata2['rfid'].isin(self.exclude)]
        excl_metadata3 = trait_metadata3[~trait_metadata3['rfid'].isin(self.exclude)]
        self.trait_md_excl = [excl_metadata1, excl_metadata2, excl_metadata3]


    # test that excluded RFIDs are successfully dropped
    def test_drop_excluded_rfids(self):
        
        for i, req in enumerate(self.full_req_list):
            self.assertEqual(list(req.trait_metadata['rfid']), 
            list(self.trait_md[i]['rfid']))
        
        for i, req in enumerate(self.excl_req_list):
            self.assertEqual(list(req.trait_metadata['rfid']), 
            list(self.trait_md_excl[i]['rfid']))


    # test that metadata & predictions successfully merge
    def test_metadata_merge(self):
        
        for i, req in enumerate(self.full_req_list):
            pd.testing.assert_frame_equal(req.trait_metadata, self.trait_md[i])

        for i, req in enumerate(self.excl_req_list):
            pd.testing.assert_frame_equal(req.trait_metadata, 
            self.trait_md_excl[i])


    # test that data are sorted by the requested trait
    def test_sorted_trait(self):
        
        for req in self.all_requests:
            self.assertEqual(
                req.trait_metadata[req.trait].iloc[0],
                max(req.trait_metadata[req.trait]))
            self.assertEqual(
                req.trait_metadata[req.trait].iloc[-1],
                min(req.trait_metadata[req.trait]))

#     # # test that available rats are broken by sex and maintained in order of prediction
# #     # def test_sex_order(self):
# #     #     pass

    # test that proposal() returns the most extreme possible rats
    def test_proposal(self):

        # save exclude list as a dictionary of rats & predictions
        unavail_dicts = []
        for req in self.full_req_list:
            trait = req.trait
            excl_df = req.trait_metadata[req.trait_metadata['rfid'].\
                isin(exclude_rfids)]
            excl_dict = excl_df.set_index('rfid')\
                [['sex', trait]].to_dict(orient='index')
            unavail_rats = \
            {k: (v['sex'], v[trait]) for k, v in excl_dict.items()}
            unavail_dicts.append(unavail_rats)

        # test that a proposal made without an exclude list
        # returns the RFIDs with the most extreme trait values
        for req in self.full_req_list:
        
            max_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
            min_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
            max_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), req.trait]
            min_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), req.trait]

            proposed_delta, proposed_rfids = req.proposal()
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

        # test that a proposal made with an exclude list
        # returns the RFIDs with the most extreme trait values
        for req in self.excl_req_list:
            
            max_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
            min_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
            max_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), req.trait]
            min_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), req.trait]

            proposed_delta, proposed_rfids = req.proposal()

            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

        # test that subsequent proposals 
        # return the RFIDs with the most extreme trait values
        for req in self.all_requests:
            
            max_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
            min_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
            max_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), req.trait]
            min_pred = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), req.trait]

            proposed_delta, proposed_rfids = req.proposal()
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)
        
            proposed_delta, proposed_rfids = req.proposal()
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

        # test explicit input of the exclude list
        for i, req in enumerate(self.full_req_list):
            
            max_rat = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nlargest(1, req.trait)['rfid'].values[0]
            min_rat = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nsmallest(1, req.trait)['rfid'].values[0]
            max_pred = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nlargest(1, req.trait)[req.trait].values[0]
            min_pred = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nsmallest(1, req.trait)[req.trait].values[0]
            
            # test using a list as input
            proposed_delta, proposed_rfids = req.proposal(exclude_rfids)
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

            # test using a dictionary as input
            proposed_delta, proposed_rfids = req.proposal(unavail_dicts[i])
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

        for i, req in enumerate(self.excl_req_list):
            
            max_rat = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nlargest(1, req.trait)['rfid'].values[0]
            min_rat = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nsmallest(1, req.trait)['rfid'].values[0]
            max_pred = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nlargest(1, req.trait)[req.trait].values[0]
            min_pred = req.trait_metadata.\
                loc[~req.trait_metadata['rfid'].isin(list(unavail_dicts[i].keys())), \
                    ['rfid', req.trait]].nsmallest(1, req.trait)[req.trait].values[0]
            
            # test using a list as input
            proposed_delta, proposed_rfids = req.proposal(exclude_rfids)
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)

            # test using a dictionary as input
            proposed_delta, proposed_rfids = req.proposal(unavail_dicts[i])
            self.assertEqual(min_rat, proposed_rfids[0])
            self.assertEqual(max_rat, proposed_rfids[1])
            self.assertEqual(max_pred - min_pred, proposed_delta)


    # test that rats are removed from the available list after assigned
    def test_remove(self):
        
        for req in self.all_requests:
                    
            # get RFIDs and sex of rats with the min & max trait values, 
            # plus one random rat
            max_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
            min_rat = req.trait_metadata.\
                loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
            random_rat = req.trait_metadata['rfid'].sample(1).values[0]

            max_rat_sex = req.trait_metadata\
                [req.trait_metadata['rfid']==max_rat]['sex'].values[0]
            min_rat_sex = req.trait_metadata\
                [req.trait_metadata['rfid']==min_rat]['sex'].values[0]
            random_rat_sex = req.trait_metadata\
                [req.trait_metadata['rfid']==random_rat]['sex'].values[0]
            
            sampled_rats = [max_rat, min_rat, random_rat]
            sampled_sexes = [max_rat_sex, min_rat_sex, random_rat_sex]
            sampled_rats = dict(zip(sampled_rats, sampled_sexes))
            for rfid, sex in sampled_rats.items():
                if sex == 'M':
                    available_by_sex = req.available_males
                elif sex == 'F':
                    available_by_sex = req.available_females

                # first ensure RFIDs are present in the list of available RFIDs
                # and in data dictionaries for available rats
                self.assertIn(rfid, req.available_rfids)
                self.assertIn(rfid, req.available_rats)
                self.assertIn(rfid, available_by_sex)
                
                # then remove sampled RFIDs and re-fetch available_by_sex
                req.remove([rfid])

                # update available_by_sex after removal
                if sex == 'M':
                    available_by_sex = req.available_males
                elif sex == 'F':
                    available_by_sex = req.available_females

                # verify RFIDs were removed from the list of available RFIDs
                # and data dictionaries for available rats
                self.assertNotIn(rfid, req.available_rfids)
                self.assertNotIn(rfid, req.available_rats)
                self.assertNotIn(rfid, available_by_sex)


    # test that RATTACA rats are properly assigned
    def test_assign_rattaca(self):

        for req in self.all_requests:
    
            # get RFIDs and sex of rats with the min & max trait values 
            md = req.trait_metadata
            md = md[md['rfid'].isin(req.available_rfids)]
            max_rat = md.loc[md[req.trait].idxmax(), 'rfid']
            min_rat = md.loc[md[req.trait].idxmin(), 'rfid']
            max_rat_sex = md[md['rfid']==max_rat]['sex'].values[0]
            min_rat_sex = md[md['rfid']==min_rat]['sex'].values[0]
            
            # manually calculate the delta that would result from assignment
            max_rat_pred = req.available_rats[max_rat][1]
            min_rat_pred = req.available_rats[min_rat][1]
            assignment_delta = max_rat_pred - min_rat_pred
            assignment_delta = req.delta + assignment_delta
            
            # assign sampled rats
            req.assign_rattaca([min_rat, max_rat])

            # check that assign_rattaca() assigns the pair that maximizes delta
            self.assertEqual(assignment_delta, req.delta)


    # test that is_satisfied_rattaca() correctly distinguishes between requests
    # that have or have not yet been filled
    def test_is_satisfied_rattaca(self):

        for req in self.all_requests:
            
            # skp project 2 (no males to assign)
            if req.n_requested_males == 0:
                continue

            # fully assign to males_high
            md = req.trait_metadata
            md_males_high = md[(md['rfid'].isin(req.available_rfids)) & \
                (md[f'{req.trait}_2group'] == 'high') & (md['sex'] == 'M')]
            males_high = md_males_high['rfid']\
                .tolist()[0 : (req.n_requested_males_high + 1)]
            extra_male = males_high.pop()

            # fully assign to females_low
            md_females_low = md[(md['rfid'].isin(req.available_rfids)) & \
                (md[f'{req.trait}_2group'] == 'low') & (md['sex'] == 'F')]
            females_low = md_females_low['rfid']\
                .tolist()[0 : (req.n_requested_females_low + 1)]
            extra_female = females_low.pop()

            self.assertIn(extra_male, req.available_rfids)
            self.assertIn(extra_female, req.available_rfids)

            # first ensure that extras are available for assignment
            self.assertIn(extra_male, req.available_rfids)
            self.assertIn(extra_female, req.available_rfids)

            # assign all low females, high males
            for i, rfid in enumerate(males_high):
                assign_pair = [females_low[i], rfid]
                print(f'assign_pair: {assign_pair}')
                print(assign_pair[0])
                req.assign_rattaca(assign_pair)
                print(f'md shape: {req.trait_metadata.shape}')
                print(f'avail rats: {len(req.available_rats)}')
                print(extra_male in req.available_rfids)

            mh_assigned = req.is_satisfied('M','high')
            ml_assigned = req.is_satisfied('M','low')
            fh_assigned = req.is_satisfied('F','high')
            fl_assigned = req.is_satisfied('F','low')

            # test complete vs incomplete assignments
            self.assertEqual(mh_assigned, True)
            self.assertEqual(ml_assigned, False)
            self.assertEqual(fh_assigned, False)
            self.assertEqual(fl_assigned, True)

            # test that extras are removed from eligibility for assignment
            self.assertNotIn(extra_male, req.available_rfids)
            self.assertNotIn(extra_female, req.available_rfids)

            # test that extras cannot be assigned
            extra_pair = [extra_female, extra_male]
            extra_assign = req.assign_rattaca(extra_pair)
            self.assertIsNone(extra_assign)

    # test that counter properties correctly track assigned samples
    def test_properties(self):

        # test that, before any assignments, all assignment counts are zero
        for req in self.all_requests:
            for sample in req.n_assigned:
                self.assertEqual(req.n_assigned[sample], 0)
            self.assertEqual(req.n_remaining['n_remaining_total'], \
                req.n_requested_total)
            self.assertEqual(len(req.available_rats),req.trait_metadata.shape[0])

    #     ### TO DO: test for updates to properties after assignments


if __name__ == '__main__':
    unittest.main()