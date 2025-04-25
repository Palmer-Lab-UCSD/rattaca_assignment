"""
Basic utility functions that don't depend on request model classes.
"""

import pandas as pd

def prep_colony_df(args):
    '''
    Drops unusable samples from the colony dataframe.
    
    Drops all dead rats, rats with undetermined sex, rats lacking RFIDs,
    and excluded RFIDs from consideration for assignment. Adds one column
    denoting which rats have been genotyped (for use in assigning to 
    random-choice projects).

    Args:
        args.colony_dataframe: The path to the colony dataframe (csv).
        args.predictions: The path to the predictions csv.
        args.exclude: The path to the exclude list (csv).

    Returns:
        A pandas dataframe with colony data for all rats available 
        for assignment, with an added columns denoting which samples 
        have been genotyped.
    '''

    df = pd.read_csv(args.colony_dataframe[0], 
                                dtype = {'rfid': str, 'accessid': int})
    
    if args.predictions:
        preds_df = pd.read_csv(args.predictions[0], dtype = {'rfid': str})
        gtyped_rfids = preds_df['rfid'].tolist()
    else:
        gtyped_rfids = []

    # check for duplicate RFIDs
    if sum(df.duplicated(subset = ['rfid'])) > 0:
        print('Error: Colony metadata has duplicate RFIDs')
        print('Check the following samples before proceeding:')
        dups = df.duplicated(subset = ['rfid'])
        dups = df[dups]
        for rfid in dups['rfid']:
            print(rfid)
        exit()

    # drop rats without RFIDs
    df[~df['rfid'].isnull()]
    
    # identify dead rats
    dead_strs = ['dead', 'die', 'death', 'euth', 'eauth', 'kill', 'starve']
    dead_search = '|'.join(dead_strs)
    dead_rats = df[df['comments']\
        .str.contains(dead_search, case=False, na=False)]['rfid'].tolist()
    
    # identify rats with unknown sex
    ambig_sex = df[df['comments']\
        .str.contains('unsure sex', case=False, na=False)]['rfid'].tolist()
    
    # identify flooded cages
    flood_strs = ['flood', 'drown']
    flood_search = '|'.join(flood_strs)
    flooded_fams = df[df['comments']\
        .str.contains(flood_search, case=False, na=False)]['breederpair'].tolist()
    flooded_rats = df[df['breederpair'].isin(flooded_fams)]['rfid'].tolist()
    
    # drop dead rats, rats w/ unknown sex, rats from flooded cages
    drop_rats = dead_rats + ambig_sex + flooded_rats
    df = df[~df['rfid'].isin(drop_rats)]

    # identify rats that have been genotyped
    df['gtyped'] = df['rfid'].isin(gtyped_rfids).astype(int)

    return(df)


def plot_assignments(preds, 
                    assignments,
                    trait, 
                    assignment_col,
                    n_groups, 
                    jitter=(80, 40),  # jitter for (background, assigned) points; set None to turn off jitter
                    gen=None, 
                    random_seed=1,
                    trait_name=None,
                    outdir=None):
    """
    Produce an assignment 'S' plot: trait rank vs. trait prediction.
    
    Parameters
    ----------
    preds : str
        File path to a csv dataframe of trait predictions including a trait 
        of interest (as produced by rattaca::get_ranks_zscores()).

    assignments : str
        File path to a csv dataframe of RATTACA request assignments (as 
        produced by rattaca_assign.assign.output_assignment_results())
        
    trait : str
        The name of the trait of interest. Must be the column header for the 
        trait predictions to be plotted.

    assignment_col : str
        The column header for the column in assignments that identifies which
        assignments to be plotted.
        
    n_groups : int
        The number of group assignments to plot. 2 will assign 'high' and 
        'low' groups, 3 will assign 'high', 'med', and 'low', higher numbers 
        will return groups numbered by quantile.
    
    jitter : tuple
        (Default (80,40)) A tuple of jitter magnitudes for non-assigned and 
        assigned points, respectively. This may need tweaking for aesthetic
        plotting.
    
    gen : int
        The RATTACA generation being plotted.
    
    random_seed : int
        (Default 1) An integer value with which to set the seed for random
        jittering. This allows consistent placement of jittered points and may 
        need tweaking for aesthetic plotting.
    
    trait_name : str
        (Default None) The desired trait name to use in the figure title and 
        file name. Use this option to simplify naming when a trait has a long 
        or unintuitive variable name.
    
    outdir : str
        (Default None) The directory path in which to save the figure. If None, 
        the figure will display using plt.show().
    
    Returns
    -------
    None
        Displays figure if outdir=None. Saves a png file if outdir is not None.
    """
    
    if trait_name is None:
        trait_name = trait

    # read in files
    preds = pd.read_csv(preds, dtype = {'rfid': str})
    assigns = pd.read_csv(assignments, dtype = {'rfid': str})
    
    # define output file path if saving
    if outdir is not None:
        if jitter is not None:
            out_stem = os.path.join(outdir, f'rattaca_gen{gen}_{trait_name}_assignment_jittered.png')
        else:
            out_stem = os.path.join(outdir, f'rattaca_gen{gen}_{trait_name}_assignment_.png')

    # columns to process for plotting
    trait_rank = f'{trait}_rank'
    group_col = f'{trait}_group'
    
    # filter for assigned rows (handle different types of True values)
    assigned_df = assigns[(assigns[assignment_col] == 1) | 
                  (assigns[assignment_col] == True) | 
                  (assigns[assignment_col] == 'True')]
    assigned_rfids = assigned_df['rfid'].tolist()

    # extract predictions & group designations for the trait of interest
    trait_cols = [col for col in preds.columns if col.startswith(trait) or col == 'rfid']
    trait_preds  = preds[trait_cols].copy()
    trait_preds[group_col] = trait_groups(preds=preds, trait=trait, n_groups=n_groups)

    
    # trait dataset: predictions and assignments
    trait_df = pd.merge(assigned_df, trait_preds, how='outer', on='rfid')
                        
    # set random seed for reproducible jitter
    np.random.seed(random_seed)
    
    # get points with jitter if needed
    if jitter is not None:
        # add jitter to all points
        x_all = trait_df[trait_rank].values + np.random.uniform(-jitter[0], jitter[0], len(trait_df))
        y_all = trait_df[trait].values + np.random.uniform(-jitter[0], jitter[0], len(trait_df))
        
        # add jitter to high trait values
        high_df = trait_df[(trait_df[group_col] == 'high') & (trait_df['rfid'].isin(assigned_rfids))]
        x_trait_high = high_df[trait_rank].values + np.random.uniform(-jitter[1], jitter[1], len(high_df))
        y_trait_high = high_df[trait].values + np.random.uniform(-jitter[1], jitter[1], len(high_df))
        
        # add jitter to low trait values
        low_df = trait_df[(trait_df[group_col] == 'low') & (trait_df['rfid'].isin(assigned_rfids))]
        x_trait_low = low_df[trait_rank].values + np.random.uniform(-jitter[1], jitter[1], len(low_df))
        y_trait_low = low_df[trait].values + np.random.uniform(-jitter[1], jitter[1], len(low_df))
    else:
        # no jitter
        x_all = trait_df[trait_rank].values
        y_all = trait_df[trait].values
        
        high_df = trait_df[(trait_df[group_col] == 'high') & (trait_df['rfid'].isin(assigned_rfids))]
        x_trait_high = high_df[trait_rank].values
        y_trait_high = high_df[trait].values
        
        low_df = trait_df[(trait_df[group_col] == 'low') & (trait_df['rfid'].isin(assigned_rfids))]
        x_trait_low = low_df[trait_rank].values
        y_trait_low = low_df[trait].values
    
    # create inferno colormap 
    inferno_cmap = plt.cm.inferno
    
    # create figure and plot
    plt.figure(figsize=(7, 5))
    
    # plot background points
    plt.scatter(x_all, y_all, s=30, color='black', alpha=0.15)
    
    # plot assigned high trait values
    plt.scatter(x_trait_high, y_trait_high, s=40, color=inferno_cmap(0.7), zorder=3)
    plt.scatter(x_trait_high, y_trait_high, s=40, facecolors='none', edgecolors='black', linewidth=1.4, zorder=4)
    
    # plot assigned low trait values
    plt.scatter(x_trait_low, y_trait_low, s=40, color=inferno_cmap(0.25), zorder=3)
    plt.scatter(x_trait_low, y_trait_low, s=40, facecolors='none', edgecolors='black', linewidth=1.4, zorder=4)
    
    # set labels and title
    plt.xlabel(f"{trait}\nprediction rank", fontweight='bold')
    plt.ylabel(f"{trait}\nprediction", fontweight='bold')
    
    if gen is not None:
        plt.title(f"RATTACA gen{gen}\n{trait_name} assignments", pad=15)
    else:
        plt.title(f"{trait_name} assignments", pad=15)
    
    plt.tight_layout()
    
    # save or display
    if outdir is not None:
        plt.savefig(out_stem, dpi=300)
        plt.close()
    else:
        plt.show()
