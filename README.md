# CVS BIDS app
For information about CVS itself, see its [paper](http://www.ajnr.org/content/early/2018/09/13/ajnr.A5765).

Usage:
```
usage: run.py [-h]
              [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
              [--skullstrip] [--thresh [THRESH]] [--n4] [--cpus CPUS]
              [--bids-filter-file FILE] [--skip_bids_validator] [-v]
              bids_dir output_dir {participant,group}

CVS entrypoint script

positional arguments:
  bids_dir              The directory with the input dataset formatted
                        according to the BIDS standard.
  output_dir            The directory where the output files should be stored.
                        If you are running group level analysis this folder
                        should be prepopulated with the results of the
                        participant level analysis.
  {participant,group}   Level of the analysis that will be performed. Multiple
                        participant level analyses can be run independently
                        (in parallel) using the same output_dir.

optional arguments:
  -h, --help            show this help message and exit
  --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        The label(s) of the participant(s) that should be
                        analyzed. The label corresponds to
                        sub-<participant_label> from the BIDS spec (so it does
                        not include "sub-"). If this parameter is not provided
                        all subjects should be analyzed. Multiple participants
                        can be specified with a space separated list.
  --skullstrip          Whether to skull strip inputs
  --thresh [THRESH]     Threshold for binary segmentation mask
  --n4                  Whether to N4 correct input
  --cpus CPUS           Number of CPUs to use
  --bids-filter-file FILE
                        a JSON file describing custom BIDS input filters using
                        PyBIDS. For further details, please check out https://
                        fmriprep.readthedocs.io/en/latest/faq.html#how-do-I-
                        select-only-certain-files-to-be-input-to-fMRIPrep
  --skip_bids_validator
                        Whether or not to perform BIDS dataset validation
  -v, --version         show program's version number and exit

```