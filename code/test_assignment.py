import rattaca_assignment as ra
import unittest
import argparse
import os
import pandas as pd


# class TestParser(unittest.TestCase):
#     def test_simple(self):
#         args = ra.parse_args(['-r','first_file','second_file'])
#         print(f'request files: {args.requests}')
#         # test that the first args element equals 'first_file'
#         self.assertEqual(args.requests[0], 'first_file')

    # def test_parse_args(self):
    #     with self.assertRaises(NotImplementedError):
    #         ra.parse_args([])

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
        self.excl_request = ra.Request(trait_1, args_excl) # uses the exclude file
        self.full_request = ra.Request(trait_1, args_all) # all rats
        self.request = self.full_request

        self.metadata = pd.read_csv(metadata, dtype = {'rfid': str, 
            'accessid': int})
        self.preds = pd.read_csv(preds, dtype = {'rfid': str})
        self.exclude = pd.read_csv(exclude, header = None).astype(str)[0].values
    
    # test that metadata & predictions successfully merge
    def test_metadata_merge(self):
        expected_merged_df = pd.merge(self.metadata, self.preds, 
                                      on='rfid', how='right')
        self.assertEqual(self.full_request.rats_metadata.shape, 
            expected_merged_df.shape)

    # test that excluded RFIDs are successfully dropped
    def test_drop_excluded_rfids(self):
        expected_merged_df = pd.merge(self.metadata, self.preds, 
                                      on='rfid', how='right')
        expected_exclusion_df = expected_merged_df[~expected_merged_df['rfid']\
            .isin(self.exclude)]
        self.assertEqual(self.excl_request.rats_metadata.shape, 
            expected_exclusion_df.shape)

    # test that data are sorted by the requested trait
    def test_sorted_trait(self):
        self.assertEqual(
            self.request.rats_metadata[self.request.trait][0],
            max(self.request.rats_metadata[self.request.trait]))

    # test that proposal() returns the most extreme possible rats
    def test_proposal(self):

       # rfids in the exclude list
        unavail_rfids = ['933000320946458','933000320946459','933000320946463',
            '933000320946464','933000320946465','933000320946466','933000320946467',
            '933000320946471','933000320946472','933000320946474']

        max_trait_all = max(self.request.rats_metadata[self.request.trait])
        min_trait_all = min(self.request.rats_metadata[self.request.trait])
        max_rat_all = max(self.request.available_trait_vals, 
                    key = self.request.available_trait_vals.get)
        min_rat_all = min(self.request.available_trait_vals, 
                    key = self.request.available_trait_vals.get)
        
        max_trait_excl = max(self.excl_request.rats_metadata[self.excl_request.trait])
        min_trait_excl = min(self.excl_request.rats_metadata[self.excl_request.trait])
        max_rat_excl = max(self.excl_request.available_trait_vals, 
                    key = self.excl_request.available_trait_vals.get)
        min_rat_excl = min(self.excl_request.available_trait_vals, 
                    key = self.excl_request.available_trait_vals.get)
    

        # test that a proposal made without an exclude list
        # returns the RFIDs with the most extreme trait values
        proposed_delta, proposed_rats = self.request.proposal(unavail_rfids)
        self.assertEqual(max_rat_all, proposed_rats[0])
        self.assertEqual(min_rat_all, proposed_rats[1])
        self.assertEqual(max_trait_all - min_trait_all, proposed_delta)

        # test that the above works when given an empty list
        proposed_delta, proposed_rats = self.request.proposal([])
        self.assertEqual(max_rat_all, proposed_rats[0])
        self.assertEqual(min_rat_all, proposed_rats[1])
        self.assertEqual(max_trait_all - min_trait_all, proposed_delta)

        # test that a proposal made with an exclude list
        # returns the RFIDs with the most extreme trait values
        proposed_delta, proposed_rats = self.excl_request.proposal(unavail_rfids)
        self.assertEqual(max_rat_excl, proposed_rats[0])
        self.assertEqual(min_rat_excl, proposed_rats[1])
        self.assertEqual(max_trait_excl - min_trait_excl, proposed_delta)

        # test that the above works when given an empty list
        proposed_delta, proposed_rats = self.excl_request.proposal([])
        self.assertEqual(max_rat_excl, proposed_rats[0])
        self.assertEqual(min_rat_excl, proposed_rats[1])
        self.assertEqual(max_trait_excl - min_trait_excl, proposed_delta)
        
    
    # # test that rats are removed from the available list after assigned
    # def test_remove(self):
    #     with self.assertRaises(NotImplementedError):
    #         self.request.remove()


    # def test_update(self):
    #     with self.assertRaises(NotImplementedError):
    #         self.request.update()



#     def test_is_satistfied(self):
#         with self.assertRaises(NotImplementedError):
#             self.request.is_satisfied()

#     def test_delta(self):
#         with self.assertRaises(NotImplementedError):
#             self.request.delta

#     # tidy up the test fixture after testing
#     def tearDown(self):
#         self.request.dispose()



if __name__ == '__main__':
    unittest.main()