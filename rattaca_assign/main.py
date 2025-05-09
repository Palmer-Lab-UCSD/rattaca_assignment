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
        
        assignments_df = output_assignment_results(
            args = args, 
            assigned_requests = assignments, 
            output_prefix = os.path.join(args.output_dir[0], args.output_prefix[0])
        )
        
        output_assignment_preds(
            assignments = assignments_df,
            preds = args.predictions[0],
            outdir = args.output_dir[0],
            args = args
        )

    else:
        print('\nNo results saved to file. Set [-o OUTPUT_PREFIX] in the command line to save results \n')

    print('\nRATTACA assignment complete! \n')


if __name__ == '__main__':
    main()