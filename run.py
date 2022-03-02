#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from bids import BIDSLayout
from glob import glob
from pathlib import Path

__version__ = open(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "version")
).read()


def collect_data(bids_dir, participant_label, filters=None, bids_validate=True):
    """
    Uses pybids to retrieve the input data for a given participant
    """
    if isinstance(bids_dir, BIDSLayout):
        layout = bids_dir
    else:
        layout = BIDSLayout(str(bids_dir), validate=bids_validate)

    queries = {
        "flair": {"datatype": "anat", "suffix": "FLAIR"},
        "t1w": {"datatype": "anat", "suffix": "T1w"},
        "epi": {"datatype": "anat", "suffix": "T2star"},
    }
    bids_filters = filters or {}
    for acq, entities in bids_filters.items():
        queries[acq].update(entities)

    subj_data = {
        dtype: sorted(
            layout.get(
                return_type="file",
                subject=participant_label,
                extension=["nii", "nii.gz"],
                **query,
            )
        )
        for dtype, query in queries.items()
    }

    return subj_data, layout


def _bids_filter(value):
    from json import loads
    from bids.layout import Query

    if value and Path(value).exists():
        try:
            filters = loads(
                Path(value).read_text(), object_hook=_filter_pybids_none_any
            )
        except Exception as e:
            raise Exception(
                "Unable to parse BIDS filter file. Check that it is " "valid JSON."
            )
    else:
        raise Exception("Unable to load BIDS filter file " + value)

    # unserialize pybids Query enum values
    for acq, _filters in filters.items():
        filters[acq] = {
            k: getattr(Query, v[7:-4])
            if not isinstance(v, Query) and "Query" in v
            else v
            for k, v in _filters.items()
        }
    return filters


def _filter_pybids_none_any(dct):
    import bids

    return {
        k: bids.layout.Query.NONE
        if v is None
        else (bids.layout.Query.ANY if v == "*" else v)
        for k, v in dct.items()
    }


def extract_run(s):
    if "run-" not in s:
        return 0
    return int(re.search(r"run-([0-9]+)", s)[1])


def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        env=merged_env,
    )
    while True:
        line = process.stdout.readline()
        line = str(line, "utf-8")[:-1]
        print(line)
        if line == "" and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d" % process.returncode)


parser = argparse.ArgumentParser(description="CVS entrypoint script")
parser.add_argument(
    "bids_dir",
    help="The directory with the input dataset "
    "formatted according to the BIDS standard.",
)
parser.add_argument(
    "output_dir",
    help="The directory where the output files "
    "should be stored. If you are running group level analysis "
    "this folder should be prepopulated with the results of the "
    "participant level analysis.",
)
parser.add_argument(
    "analysis_level",
    help="Level of the analysis that will be performed. "
    "Multiple participant level analyses can be run independently "
    "(in parallel) using the same output_dir.",
    choices=["participant", "group"],
)
parser.add_argument(
    "--participant_label",
    help="The label(s) of the participant(s) that should be analyzed. The label "
    "corresponds to sub-<participant_label> from the BIDS spec "
    '(so it does not include "sub-"). If this parameter is not '
    "provided all subjects should be analyzed. Multiple "
    "participants can be specified with a space separated list.",
    nargs="+",
)
parser.add_argument(
    "--skullstrip",
    help="Whether to skull strip inputs ",
    action="store_true",
)
parser.add_argument(
    "--thresh", help="Threshold for binary segmentation mask", nargs="?", default=0.2
)
parser.add_argument("--n4", help="Whether to N4 correct input", action="store_true")
parser.add_argument("--cpus", help="Number of CPUs to use", type=str, default="1")
parser.add_argument(
    "--bids-filter-file",
    dest="bids_filters",
    action="store",
    type=_bids_filter,
    metavar="FILE",
    help="a JSON file describing custom BIDS input filters using PyBIDS. "
    "For further details, please check out "
    "https://fmriprep.readthedocs.io/en/latest/faq.html#"
    "how-do-I-select-only-certain-files-to-be-input-to-fMRIPrep",
)
parser.add_argument(
    "--skip_bids_validator",
    help="Whether or not to perform BIDS dataset validation",
    action="store_true",
)
parser.add_argument(
    "-v", "--version", action="version", version="CVS version {}".format(__version__)
)


args = parser.parse_args()

if not args.skip_bids_validator:
    run("bids-validator %s" % args.bids_dir)

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

# running participant level
if args.analysis_level == "participant":

    for subject_label in subjects_to_analyze:
        subject_data, layout = collect_data(
            args.bids_dir,
            subject_label
            if not subject_label.startswith("sub-")
            else subject_label[4:],
            filters=args.bids_filters,
            bids_validate=not args.skip_bids_validator,
        )
        t1 = sorted(subject_data["t1w"])
        flair = sorted(subject_data["flair"])
        epi = sorted(subject_data["epi"])

        if len(t1) == 0:
            raise Exception("No T1w images found for participant %s" % subject_label)
        if len(flair) == 0:
            raise Exception("No FLAIR images found for participant %s" % subject_label)
        if len(epi) == 0:
            raise Exception("No T2star images found for participant %s" % subject_label)

        sorted_t1s = []
        sorted_flairs = []
        sorted_epis = []
        layout_sessions = layout.get_sessions()
        if len(layout_sessions) == 0:
            # sort in reverse order because the last run is usually the highest quality
            t1s = sorted(t1, reverse=True, key=extract_run)
            flairs = sorted(flair, reverse=True, key=extract_run)
            epis = sorted(epi, reverse=True, key=extract_run)
            smaller = min(len(t1s), len(flairs), len(epis))
            for i in range(smaller):
                sorted_t1s.append(t1s[i])
                sorted_flairs.append(flairs[i])
                sorted_epis.append(epis[i])
                print("Pairing", t1s[i], "with", flairs[i], "and", epis[i])
        else:
            for ses in layout_sessions:
                ses = "ses-" + ses
                # sort in reverse order because the last run is usually the highest quality
                ses_t1s = sorted(
                    [t1_file for t1_file in t1 if ses in t1_file],
                    reverse=True,
                    key=extract_run,
                )
                ses_flairs = sorted(
                    [flair_file for flair_file in flair if ses in flair_file],
                    reverse=True,
                    key=extract_run,
                )
                ses_epis = sorted(
                    [epi_file for epi_file in epi if ses in epi_file],
                    reverse=True,
                    key=extract_run,
                )
                smaller = min(len(ses_t1s), len(ses_flairs), len(ses_epis))
                for i in range(smaller):
                    sorted_t1s.append(ses_t1s[i])
                    sorted_flairs.append(ses_flairs[i])
                    sorted_epis.append(ses_epis[i])
                    print(
                        "Pairing", ses_t1s[i], "with", ses_flairs[i], "and", ses_epis[i]
                    )
            if len(sorted_t1s) == 0:
                raise Exception(
                    "Subject %s has at least 1 T1w and FLAIR, but not in the same session"
                    % subject_label
                )

        outdir_suffix = ""
        for i in range(len(sorted_t1s)):
            if len(sorted_t1s) > 1:
                if "ses-" in sorted_t1s[i]:
                    outdir_suffix = (
                        "_" + re.search(r"ses-[0-9a-zA-Z]+", sorted_t1s[i])[0]
                    )
                    if "run-" in sorted_t1s[i]:
                        outdir_suffix += (
                            "_" + re.search(r"run-[0-9]+", sorted_t1s[i])[0]
                        )
                else:
                    outdir_suffix = "_%d" % i
            out_file = os.path.split(sorted_t1s[i])[-1].replace("_T1w.", "_cvs.")
            Path(args.output_dir + outdir_suffix).mkdir(parents=True, exist_ok=True)
            cmd = "/process_subject.R --outdir %s --indir %s --flair %s --t1 %s --epi %s --thresh %s" % (
                os.path.realpath(args.output_dir) + outdir_suffix,
                os.path.dirname(sorted_t1s[i]),
                os.path.basename(sorted_flairs[i]),
                os.path.basename(sorted_t1s[i]),
                os.path.basename(sorted_epis[i]),
                args.thresh,
            )
            if args.n4:
                cmd += " --n4"
            if args.skullstrip:
                cmd += " --skullstrip"
            print(cmd)
            os.chdir("/") # R script sources things using relative path
            run(
                cmd,
                {
                    "CORES": args.cpus,  # the variable the program actually uses
                    # what various libraries use
                    "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS": args.cpus,
                    "OMP_NUM_THREADS": args.cpus,
                    "OMP_THREAD_LIMIT": args.cpus,
                    "MKL_NUM_THREADS": args.cpus,
                    "OPENBLAS_NUM_THREADS": args.cpus,
                },
            )
else:
    print("CVS only supports --analysis_level=participant")
