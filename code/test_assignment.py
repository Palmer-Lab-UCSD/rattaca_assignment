import rattaca_assignment as ra
import unittest
import argparse
import os
import pandas as pd
import numpy as np
import random

### FINISH: test that group assignment works
def test_trait_groups(df, trait, n_groups=2):

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

# # tests for argument parsing
# class TestParser(unittest.TestCase):
    
#     # test that args elements are read in correctly
#     def test_parse_args_simple(self):
#         args = ra.parse_args(['-c', 'colony_df_file',
#                               '-p', 'preds_file',
#                               '-r','first_request_file','second_request_file'])
#         self.assertEqual(args.colony_dataframe[0], 'colony_df_file')
#         self.assertEqual(args.predictions[0], 'preds_file')
#         self.assertEqual(args.requests[0], 'first_request_file')
#         self.assertEqual(args.requests[1], 'second_request_file')    
    
    # # FINISH: more detailed argument parsing with real files
    # def test_parse_args(self):
    #     with self.assertRaises(NotImplementedError):
    #         ra.parse_args(['-c', 'colony_df_file',
    #                        '-p', 'preds_file'])


# tests for the Request class
class TestRequest(unittest.TestCase):

    # set up the test fixture
    def setUp(self):

        # test dataset
        proj_dir = os.path.join('../') 
        req_dir = os.path.join(proj_dir, 'proj_requests')
        data_dir = os.path.join(proj_dir, 'data')
        metadata = os.path.join(data_dir, 'hsw_gen100_colony.csv')
        preds = os.path.join(data_dir, 'rattaca_gen100_preds.csv')
        exclude = os.path.join(data_dir, 'rattaca_gen100_exclude.csv')
        trait_1 = os.path.join(req_dir, 'rattaca', 'trait_1.json')
        trait_2 = os.path.join(req_dir, 'rattaca', 'trait_2.json')
        trait_3 = os.path.join(req_dir, 'rattaca', 'trait_3.json')

        # input arguments
        input_args = ['-c', metadata, '-p', preds, '-r', trait_1, trait_2, 
                      trait_3, '-e', exclude]
        args_excl = ra.parse_args(input_args) # uses the list of RFIDs to exclude
        args_all = ra.parse_args(input_args[:-2]) # all rats

        # create Request instances for testing
        self.request1 = ra.Request(trait_1, args_all) # all rats
        self.request2 = ra.Request(trait_2, args_all)
        self.request3 = ra.Request(trait_3, args_all)
        self.excl_req1 = ra.Request(trait_1, args_excl) # uses the exclude file
        self.excl_req2 = ra.Request(trait_2, args_excl) # uses the exclude file
        self.excl_req3 = ra.Request(trait_3, args_excl) # uses the exclude file

        # read in data from arguments
        self.colony_df = pd.read_csv(metadata, dtype = {'rfid': str, 
            'accessid': int})
        self.preds_df = pd.read_csv(preds, dtype = {'rfid': str})
        self.exclude = pd.read_csv(exclude, header = None).astype(str)[0].values

        # set up metadata for each request
        self.trait_metadata1 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_1']], on='rfid', how='right')\
            .sort_values(by='trait_1', axis=0, ascending=False)\
            .dropna(subset='trait_1', ignore_index=True)
        self.trait_metadata2 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_2']], on='rfid', how='right')\
            .sort_values(by='trait_2', axis=0, ascending=False)\
            .dropna(subset='trait_2', ignore_index=True)
        self.trait_metadata3 = self.colony_df.\
            merge(self.preds_df[['rfid', 'trait_3']], on='rfid', how='right')\
            .sort_values(by='trait_3', axis=0, ascending=False)\
            .dropna(subset='trait_3', ignore_index=True)
        
        # add group assignments to each request's metadata
        self.trait_metadata1['trait_1_2group'] = \
            test_trait_groups(self.trait_metadata1, 'trait_1', n_groups=2)
        self.trait_metadata2['trait_2_2group'] = \
            test_trait_groups(self.trait_metadata2, 'trait_2', n_groups=2)
        self.trait_metadata3['trait_3_2group'] = \
            test_trait_groups(self.trait_metadata3, 'trait_3', n_groups=2)

    
#     # test that metadata & predictions successfully merge
#     def test_metadata_merge(self):
#         expected_merged_df = \
#             pd.merge(self.colony_df, self.preds_df[['rfid','trait_1']], 
#                      on='rfid', how='right')
#         # add a column of high/low group assignments
#         expected_merged_df[f'{self.request1.trait}_2group'] = \
#             test_trait_groups(self.trait_metadata1, self.request1.trait, n_groups=2)

#         self.assertEqual(list(self.request1.trait_metadata.columns), 
#             list(expected_merged_df.columns))
#         self.assertEqual(self.request1.trait_metadata.shape, 
#             expected_merged_df.shape)
        

#     # test that excluded RFIDs are successfully dropped
#     def test_drop_excluded_rfids(self):
#         expected_merged_df = \
#             pd.merge(self.colony_df, self.preds_df[['rfid','trait_1']], 
#                      on='rfid', how='right')
#         # add a column of high/low group assignments
#         expected_merged_df[f'{self.request1.trait}_2group'] = \
#             test_trait_groups(self.trait_metadata1, 
#                               self.request1.trait, n_groups=2)
#         # sort by the requested trait, remove NAs
#         expected_merged_df = expected_merged_df.sort_values(by = 'trait_1', 
#             axis=0, ascending=False)\
#                 .dropna(subset='trait_1', ignore_index=True)
#         # drop any RFIDs from the exclusion list
#         expected_merged_df = \
#             expected_merged_df\
#                 [~expected_merged_df['rfid'].isin(self.exclude)]

#         expected_available_rfids = expected_merged_df['rfid'].tolist()
#         expected_available_rfids = [rfid for rfid in expected_available_rfids
#                                 if rfid not in self.exclude]
#         self.assertEqual(expected_available_rfids, 
#                          self.excl_req1.available_rfids)


#     # test that data are sorted by the requested trait
#     def test_sorted_trait(self):

#         full_req_list = [self.request1, self.request2, self.request3]
#         excl_req_list = [self.excl_req1, self.excl_req2, self.excl_req3]
#         all_requests = full_req_list + excl_req_list
        
#         for req in all_requests:
#             self.assertEqual(
#                 req.trait_metadata[req.trait].iloc[0],
#                 max(req.trait_metadata[req.trait]))
#             self.assertEqual(
#                 req.trait_metadata[req.trait].iloc[-1],
#                 min(req.trait_metadata[req.trait]))

#     # # test that available rats are broken by sex and maintained in order of prediction
# #     # def test_sex_order(self):
# #     #     pass

#     # test that proposal() returns the most extreme possible rats
#     def test_proposal(self):

#        # rfids in the exclude list
#         unavail_rats1 = {
#             '933000320947017':('M',-0.29767011),'933000320978275':('M',0.34185398),
#             '933000320978314':('F',0.08600022),'933000320978345':('M',-0.06004481),
#             '933000320980698':('F',0.10041962),'933000320946466':('F',-0.14511259,),
#             '933000320946467':('F',-0.16900296),'933000320978230':('M',0.06880568)}
#         unavail_rats2 = {
#             '933000320947017':('M',-0.03467341),'933000320978275':('M',-0.20838943),
#             '933000320978314':('F',-0.47193199),'933000320978345':('M',0.59633627),
#             '933000320980698':('F',-0.34553970),'933000320946466':('F',0.23791499),
#             '933000320946467':('F',0.09621416),'933000320978230':('M',-0.38332357)}
#         unavail_rats3 = {
#             '933000320947017':('M',-0.0571968),'933000320978275':('M',0.04379969),
#             '933000320978314':('F',-0.00011618),'933000320978345':('M',-0.01265264),
#             '933000320980698':('F',-0.15038199),'933000320946466':('F',-0.10379708),
#             '933000320946467':('F',-0.01009846),'933000320978230':('M',0.07488263)}

#         full_req_list = [self.request1, self.request2, self.request3]
#         excl_req_list = [self.excl_req1, self.excl_req2, self.excl_req3]
#         all_requests = full_req_list + excl_req_list
#         unavail_dicts = [unavail_rats1, unavail_rats2, unavail_rats3]
#         unavail_list = list(unavail_rats1.keys())
        
#         # test that a proposal made without an exclude list
#         # returns the RFIDs with the most extreme trait values
#         for req in full_req_list:
        
#             max_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
#             min_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
#             max_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), req.trait]
#             min_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), req.trait]

#             proposed_delta, proposed_rfids = req.proposal()
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#         # test that a proposal made with an exclude list
#         # returns the RFIDs with the most extreme trait values
#         for i, req in enumerate(excl_req_list):
            
#             max_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
#             min_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
#             max_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), req.trait]
#             min_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), req.trait]

#             proposed_delta, proposed_rfids = req.proposal(unavail_dicts[i])
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#             proposed_delta, proposed_rfids = req.proposal(unavail_dicts[i])
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#         # test that subsequent proposals 
#         # return the RFIDs with the most extreme trait values
#         for req in all_requests:
            
#             max_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
#             min_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
#             max_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), req.trait]
#             min_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), req.trait]

#             proposed_delta, proposed_rfids = req.proposal()
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)
        
#             proposed_delta, proposed_rfids = req.proposal()
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#         # test using a list as input
#         for req in all_requests:
            
#             max_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
#             min_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
#             max_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), req.trait]
#             min_pred = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), req.trait]
            
#             proposed_delta, proposed_rfids = req.proposal()
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#         # test that the max/min preds are returned when excluding RFIDs
#         for req in all_requests:
            
#             max_rat = req.trait_metadata.\
#                 loc[~req.trait_metadata['rfid'].isin(unavail_list), \
#                     ['rfid', req.trait]].nlargest(1, req.trait)['rfid'].values[0]
#             min_rat = req.trait_metadata.\
#                 loc[~req.trait_metadata['rfid'].isin(unavail_list), \
#                     ['rfid', req.trait]].nsmallest(1, req.trait)['rfid'].values[0]
#             max_pred = req.trait_metadata.\
#                 loc[~req.trait_metadata['rfid'].isin(unavail_list), \
#                     ['rfid', req.trait]].nlargest(1, req.trait)[req.trait].values[0]
#             min_pred = req.trait_metadata.\
#                 loc[~req.trait_metadata['rfid'].isin(unavail_list), \
#                     ['rfid', req.trait]].nsmallest(1, req.trait)[req.trait].values[0]
            
#             proposed_delta, proposed_rfids = req.proposal(unavail_list)
#             self.assertEqual(min_rat, proposed_rfids[0])
#             self.assertEqual(max_rat, proposed_rfids[1])
#             self.assertEqual(max_pred - min_pred, proposed_delta)

#     # test that rats are removed from the available list after assigned
#     def test_remove(self):
        
#         request_list = [self.request1, self.request2, self.request3,
#                          self.excl_req1, self.excl_req2, self.excl_req3]

#         for req in request_list:
                    
#             # get RFIDs and sex of rats with the min & max trait values, 
#             # plus one random rat
#             max_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmax(), 'rfid']
#             min_rat = req.trait_metadata.\
#                 loc[req.trait_metadata[req.trait].idxmin(), 'rfid']
#             random_rat = random.randrange(0,len(req.available_rfids))
#             random_rat = req.trait_metadata['rfid'][random_rat]
#             max_rat_sex = req.trait_metadata\
#                 [req.trait_metadata['rfid']==max_rat]['sex'].values[0]
#             min_rat_sex = req.trait_metadata\
#                 [req.trait_metadata['rfid']==min_rat]['sex'].values[0]
#             random_rat_sex = req.trait_metadata\
#                 [req.trait_metadata['rfid']==random_rat]['sex'].values[0]
            
#             sampled_rats = [max_rat, min_rat, random_rat]
#             sampled_sexes = [max_rat_sex, min_rat_sex, random_rat_sex]
#             sampled_rats = dict(zip(sampled_rats, sampled_sexes))
#             for rfid, sex in sampled_rats.items():
#                 if sex == 'M':
#                     available_by_sex = req.available_males
#                 elif sex == 'F':
#                     available_by_sex = req.available_females

#                 # first ensure RFIDs are present in the list of available RFIDs
#                 # and in data dictionaries for available rats
#                 self.assertIn(rfid, req.available_rfids)
#                 self.assertIn(rfid, req.available_rats)
#                 self.assertIn(rfid, available_by_sex)
                
#                 # then remove sampled RFIDs and re-fetch available_by_sex
#                 req.remove([rfid])

#                 # update available_by_sex after removal
#                 if sex == 'M':
#                     available_by_sex = req.available_males
#                 elif sex == 'F':
#                     available_by_sex = req.available_females

#                 # verify RFIDs were removed from the list of available RFIDs
#                 # and data dictionaries for available rats
#                 self.assertNotIn(rfid, req.available_rfids)
#                 self.assertNotIn(rfid, req.available_rats)
#                 self.assertNotIn(rfid, available_by_sex)

    # test that RATTACA rats are properly assigned
    def test_assign_rattaca(self):

        request_list = [self.request1, self.request2, self.request3, 
                        self.excl_req1, self.excl_req2, self.excl_req3]

        for req in request_list:
    
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

    # def test_properties(self):

    #     request_list = [self.request1, self.request2, self.request3, 
    #             self.excl_req1, self.excl_req2, self.excl_req3]

    #     # test that, before any assignments, all assignment counts are zero
    #     for req in request_list:
    #         for key in req.n_assigned:
    #             self.assertEqual(req.n_assigned[key], 0)
    #         self.assertEqual(req.n_remaining['n_remaining_total'], \
    #             req.n_requested_total)
    #         self.assertEqual(len(req.available_rats),req.trait_metadata.shape[0])

    #     ### TO DO: test for updates to properties after assignments

    # test that is_satisfied_rattaca() correctly distinguishes between requests
    # that have or have not yet been filled
    def test_is_satisfied_rattaca(self):

        request_list = [self.request1, self.request3, 
                        self.excl_req1, self.excl_req3]

        for req in request_list:
            
            # fully assign to males_high
            md = req.trait_metadata
            md_males_high = md[(md['rfid'].isin(req.available_rfids)) & \
                (md[f'{req.trait}_2group'] == 'high') & (md['sex'] == 'M')]
            males_high = md_males_high['rfid']\
                .tolist()[0 : (req.n_requested_males_high + 1)]
            extra_male = males_high.pop()
            md_females_low = md[(md['rfid'].isin(req.available_rfids)) & \
                (md[f'{req.trait}_2group'] == 'low') & (md['sex'] == 'F')]
            females_low = md_females_low['rfid']\
                .tolist()[0 : (req.n_requested_males_high + 1)]
            extra_female = females_low.pop()

            for i, rfid in enumerate(males_high):
                assign_pair = [females_low[i], rfid]
                req.assign_rattaca(assign_pair)

            mh_assigned = req.is_satisfied_rattaca('M','high')
            ml_assigned = req.is_satisfied_rattaca('M','low')
            fh_assigned = req.is_satisfied_rattaca('F','high')
            fl_assigned = req.is_satisfied_rattaca('F','low')

            # test complete vs incomplete assignments
            self.assertEqual(mh_assigned, True)
            self.assertEqual(ml_assigned, False)
            self.assertEqual(fh_assigned, False)
            self.assertEqual(fl_assigned, True)

            # test that extras cannot be assigned
            extra_pair = [extra_female, extra_male]
            extra_assign = req.assign_rattaca(extra_pair)
            self.assertIsNone(extra_assign)


    # def test_assign(self):
    #     with self.assertRaises(NotImplementedError):
    #         self.request1.update()



#     def test_is_satistfied(self):
#         with self.assertRaises(NotImplementedError):
#             self.request1.is_satisfied()

#     def test_delta(self):
#         with self.assertRaises(NotImplementedError):
#             self.request1.delta

    # # tidy up the test fixture after testing
    # def tearDown(self):
    #     self.request1.dispose()



# if __name__ == '__main__':
#     unittest.main()