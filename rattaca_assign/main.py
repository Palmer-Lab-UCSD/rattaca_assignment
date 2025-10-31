from rattaca_assign.cli import parse_args
from rattaca_assign.core.assign import *
import os

def main(args=None):
    '''Main entry point for the RATTACA assignment package.'''

    # parse command line arguments
    args = parse_args(args)

    print('\n----------------------------')
    print('---- RATTACA Assignment ----')
    print('---------------------------- \n')

    do_step = args.step

    # if conducting new assignments
    if do_step == 'new_assignments':

        print(f'Request files:')
        for req in args.requests:
            print(f'\t{req}')
        print(f'Colony dataframe: {args.colony_dataframe[0]}')
        if args.predictions:
            print(f'Predictions file: {args.predictions[0]}')
        print('\n')
  
        # run the assignment algorithm
        assignments = run_assignments(args)

        # save results to file
        if args.output_prefix:
            
            # assignment files: per-request and all requests
            assignments_df = output_assignment_results(
                args = args, 
                assigned_requests = assignments, 
                output_prefix = os.path.join(args.output_dir[0], \
                    args.output_prefix[0])
            )
            
            # RATTACA requests results files: assignments + predictions
            output_assignment_preds(
                requests = args.requests,
                assignments = assignments_df,
                preds = args.predictions[0],
                outdir = args.output_dir[0],
                request_map = args.request_map[0],
                summary = args.preds_summary[0])


        else:
            print('\nNo results saved to file. Set [-o OUTPUT_PREFIX] in the command line to save results \n')

        print('\nRATTACA assignment complete! \n')


    # if updating previously proposed assignments
    elif do_step == 'update_assignments':

        print('----------------------------')
        print('--- Updating Assignments ---')
        print('---------------------------- \n')
    
        print(f'Assignments file: {args.assignments[0]}')
        print(f'Updates file: {args.updates[0]} \n')

        output_dir = args.output_dir[0] if args.output_dir else None
        print(f'output_dir: {output_dir}')

        update_assignments(
            assignments = args.assignments[0],
            updates = args.updates[0],
            all_requests = True,
            outdir = output_dir)
    

    # if formatting final outputs for RATTACA assignments
    elif do_step == 'rattaca_results':

        print('-----------------------------')
        print('--- Final RATTACA Results ---')
        print('----------------------------- \n')

        print(f'Assignments file: {args.assignments[0]}')

        # RATTACA requests results files: assignments + predictions
        output_assignment_preds(
            requests = args.requests,
            assignments = args.assignments[0],
            preds = args.predictions[0],
            outdir = args.output_dir[0],
            request_map = args.request_map[0],
            summary = args.preds_summary[0])



if __name__ == '__main__':
    main()