# RATTACA assignment :rat: :dart:

A Python package for assigning HS West rats to RATTACA projects based on genetic predictions.  

:construction: :construction: Under Construction :construction: :construction:

[RATTACA](https://ratgenes.org/rattaca/) is a research service that (1) uses genome-wide SNP genotypes to predict phenotypes in HS rats, then (2) uses these trait predictions to assign for study animals with predicted extreme phenotypes. This package provides helper tools for step (2): assigning HS rats to projects requesting phenotypically extreme animals. 

## How It Works

At the heart of the `rattaca_assign` package is a "constrained resource allocation" algorithm. The algorithm assigns limited resources (rats) to different projects while satisfying multiple constraints (trait prediction, sex, litter, age, etc.). RATTACA uses a greedy strategy (selecting extreme-valued animals first) and permutation to find the assignment order that maximizes a the differences in trait values assigned to all receiving projects.  

The goal of RATTACA is to provide investigators with samples of rats that are expected to statistically differ in a given trait of interest. This ideally entails assigning to a a given request all of the most extreme rats from the current cohort (two samples with high and low trait predictions). This indexing of extreme trait values for assignment is one of the core functionalities of the `rattaca_assign` package. 

In some cases, if multiple investigators request rats from a cohort based on predictions for correlated traits, conflicts over assignments must be resolved using a tiebreaking algorithm. In these cases, RATTACA assigns rats using permutation to maximize per-project, between-sample differences in predicted traits. `rattaca_assign` orders projects (decides who gets "dibs" first, second, etc), simulates assignments to each project, calculates the magnitude difference in trait scores for each simulation, then permutes alternative proposals over each possible ordering of projects. Final assignments are made in the order that maximizes the overall magnitude difference across conflicting requests. This ensures each request is allotted a maximally-extreme sample of rats (within constraints imposed by other projects) that is still expected to produce statistically distinct sub-samples.


The main assignment steps are:

1. Read request files to determine which types of assignments are requested.
2. Process colony data to identify rats available for assignment.
3. For HSW breeder requests, assign priority breeders (a single male or female from a given litter) first.
4. For RATTACA requests, try all possible permutations of request order.
5. For each permutation, propose assignments and calculate the total difference in trait values.
6. Choose the permutation that maximizes the total delta.
7. Assign rats according to the best permutation.
8. Assign breeders when only one rat (per sex) remains from a litter
9. Assign random requests only as the number of available litters <= than the number of rats requested (per sex)

Steps 8 and 9 of the algorithm ensure that extreme rats are prioritized for assignment to RATTACA requests with minimal conflict with non-RATTACA projects.

## Installation

RATTACA assignments are conducted in Python. You will need to install Python >= v3.6 with the following packages:
- pandas >=2.21 (to manipulate dataframes)
- numpy >=1.26.4 (for quantiles and binning trait predictions)

To install the package, clone the GitHub repository and install using pip:
```bash
#!/bin/bash

# clone the repository
git git@github.com:Palmer-Lab-UCSD/rattaca_assignment.git
cd rattaca_assignment.git

# install the package
pip install -e .
```

## Usage

### Command Line Interface

RATTACA assignment can be run in full from the command line:

```bash
#!/bin/bash 

rattaca_assign new_assignments -c colony_dataframe.csv -p predictions.csv -m request_map.csv -s predictions_summary.csv -o /output/path/ -f output_prefix -r request1.json request2.json
```

Command line options:
- `step`: A positional argument to define which step in the assignment process to execute. Options include
  * `new_assignments`: Assigns rats and saves formatted results.
  * `rattaca_results`: Saves formatted results for all requests given a multi-request assignment file. Useful following manual assignments.
  * `update_assignments`: Re-formats results to reflect any changes in assignments reflected in shipping sheets. Useful when executed assignments (actual rats assigned and shipped) differ from assignments previously proposed using `new_assignments`.
- `-c, --colony_df`: Path to the HS West colony dataframe CSV for the current generation (required for new assignments)
- `-r, --request_files`: Path to one or more request files in JSON format (required for new assignments)
- `-p, --predictions`: Path to the RATTACA predictions CSV (optional for new assignments)
- `-o, --output_dir`: Output directory (optional for new assignments)
- `-f, --output_prefix`: Output filename prefix (optional for new assignments)
- `-e, --exclude_rfids`: Path to CSV file containing RFIDs to exclude from assignment (optional for new assignments)
- `-a, --assignments`: Path to CSV file containing previously proposed assignments (required for updating assignments)
- `-u, --updates`: Path to JSON file specifying shipping sheets for assigned requests (required for updating assignments)
- `-s, --predictions_summary`: Path to CSV 'predictions summary' file produced during trait prediction (required for new or updated results formatting). 
- `-m, --request_map`: Path to CSV file outlining which trait predictions to include in saved results for each request (required for new and updated output formatting)
- `-v --verbose`: An optional flag to include if verbose output is desired (as for debugging).

To format per-request results following manual assignments:
```bash
#!/bin/bash 

rattaca_assign rattaca_results -a all_assignments.csv -p predictions.csv -m request_map.csv -s predictions_summary.csv -o /output/path/ -r request1.json request2.json
```
(note that no output prefix is needed. It will be taken automatically from input files)

To produce updated results following animal shipments:
```bash
#!/bin/bash 

rattaca_assign update_assignments -a all_assignments.csv -u assignment_updates.json -o /output/path/ 
```
(note that no output prefix is needed. It will be taken automatically from input files)

### Python API

The package can also be used programmatically:

```python
from rattaca_assign.core.assignment import run_assignments
from rattaca_assign.core.utils import prep_colony_df

# parse command line args
args = parse_args(['--colony_df', 'colony.csv', '--predictions', 'preds.csv', '--request_files', 'request1.json', 'request2.json'])

# run the assignment algorithm
assignments = run_assignments(args)

# output the results
for project, rfids in assignments.items():
    print(f"Project: {project}, RFIDs: {rfids}")
```

### Testing

To run unit tests for all package modules, run this bash command from the package top directory:

```bash
#!/bin/bash

python -m unittest discover
```

## Request Files

`rattaca_assign` requires at least one request file as an input. Request files are JSON files that specify a given project's request details and assignment parameters. They must be formatted depending on project type, as outlined below. Three types of request files are supported:

1. RATTACA requests - Assignments based on genetic predictions
2. Random requests - Assignments made semi-randomly, usually within some constraints (such as a limit on the number of siblings allowed from a given family)
3. HSW Breeder requests - Assignments to the breeder pool at the HS West colony.

Example json file for a RATTACA project requesting 10 male and 10 female rats with high or low predictions for the hypothetical trait "anxiety_score". The trait name corresponds to a column in `predictions.csv`. `min_age` and `max_age` are set to constrain the desired age range for the project, and `receive_date` sets the specific date (in YYYMMDD format) on which the age restrictions apply:

```json
{
  "assignment_type": "rattaca",
  "project": "AnxietyProject",
  "request_name": "rattaca_gen100_anxiety",
  "assignment_name": "AnxProj_anxiety",
  "trait": "anxiety_score",
  "min_age": 40,
  "max_age": 50,
  "receive_date": 20250101,
  "max_per_sex": 1,
  "n_rats": {
    "total": 20,
    "male": {
      "total": 10,
      "high": 5,
      "low": 5
    },
    "female": {
      "total": 10,
      "high": 5,
      "low": 5
    }
  }
}
```

Example json files for "random" request. Here, trait names are null, as assignment does not rely on any trait predictions. Often, a given project may wish to distinguish between randomly-assigned experimental and control groups, in which case two files should be made with the the same project name but different request names (and relevant request numbers per group). Here, this example will allow no same-sex siblings to be drawn from the same litter in the experimental group, but has no such restriction for the control group. Note that age restrictions and the receipt date are all null, meaning no restrictions are set:
```json
{
  "assignment_type": "random",
  "project": "AddictionProject",
  "request_name": "hsw_gen100_addiction_experimental",
  "assignment_name": "addproj_addiction_exp",
  "trait": null,
  "max_males_per_fam" : 1,
  "max_females_per_fam" : 1,
  "min_age": null,
  "max_age": null,
  "receive_date": null,
  "n_rats": {
    "total": 16,
    "male": {
      "total": 8
    },
    "female": {
      "total": 8
    }
  }
}

{
  "assignment_type": "random",
  "project": "AddictionProject",
  "request_name": "hsw_gen100_addiction_control",
  "assignment_name": "addproj_addiction_ctrl",
  "trait": null,
  "max_males_per_fam" : 0,
  "max_females_per_fam" : 0,
  "min_age": null,
  "max_age": null,
  "receive_date": null,
  "n_rats": {
    "total": 16,
    "male": {
      "total": 8
    },
    "female": {
      "total": 8
    }
  }
}
```

Example json file for HSW breeders. Assignment type, project name, request name, and assignment name are always "hsw_breeders", and trait name is always null. Typical generations at HSW keep 70 breeding pairs - one male and one female per litter:
```json
{
  "assignment_type": "hsw_breeders",
  "project": "hs_west_colony",
  "request_name": "hsw_gen100_colony",
  "assignment_name": "hsw_breeders",
  "trait": null,
  "max_males_per_fam" : 1,
  "max_females_per_fam" : 1,
  "min_age": null,
  "max_age": null,
  "receive_date": null,
  "n_rats": {
    "total": 140,
    "male": {
      "total": 70
    },
    "female": {
      "total": 70
    }
  }
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.