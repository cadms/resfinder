# python3 /home/zack/development/resfinder/run_resfinder.py -ifa 1237-49_S49.raw.fasta -acq -l 0.6 -t 0.9

# !/usr/bin/env python3
import sys
import os
from argparse import ArgumentParser
import time

start_time = time.time()

parser = ArgumentParser(allow_abbrev=False)

# General options
parser.add_argument("-ifd",
                    help="Input fasta directory.",
                    default=None)
parser.add_argument("-o",
                    help="Path to blast output",
                    default='')
parser.add_argument("-b",
                    help="Path to blastn",
                    default='blastn')
parser.add_argument("-k",
                    help="Path to KMA",
                    default=None)
parser.add_argument("-s",
                    help="Species in the sample",
                    default=None)
parser.add_argument("-l",
                    help="Minimum (breadth-of) coverage",
                    type=float,
                    default=0.60)
parser.add_argument("-t",
                    help="Threshold for identity",
                    type=float,
                    default=0.80)

# Acquired resistance options
parser.add_argument("-db_res",
                    help="Path to the databases for ResFinder",
                    default=None)
parser.add_argument("-db_res_kma",
                    help="Path to the ResFinder databases indexed with KMA. \
                          Defaults to the 'kma_indexing' directory inside the \
                          given database directory.",
                    default=None)
parser.add_argument("-d",
                    help="Databases chosen to search in - if none is specified\
                          all is used",
                    default=None)
parser.add_argument("-acq",
                    action="store_true",
                    help="Run resfinder for acquired resistance genes",
                    default=False)
parser.add_argument("-ao",
                    help="Genes are allowed to overlap this number of\
                          nucleotides. Default: 30.",
                    type=int,
                    default=30)

# Point resistance option
parser.add_argument("-c",
                    action="store_true",
                    help="Run pointfinder for chromosomal mutations",
                    default=False)
parser.add_argument("-db_point",
                    help="Path to the databases for PointFinder",
                    default=None)
parser.add_argument("-g",
                    nargs='+',
                    help="Specify genes existing in the database to \
                          search for - if none is specified all genes are \
                          included in the search.",
                    default=None)
parser.add_argument("-u",
                    action="store_true",
                    help="Show all mutations found even if in unknown to the\
                          resistance database",
                    default=False)

# # Temporary option only available temporary
# parser.add_argument("--pickle",
#                     action="store_true",
#                     help="Create a pickle dump of the Isolate object. \
#                           Currently needed in the CGE webserver. Dependency \
#                           and this option is being removed.",
#                     default=False)

args = parser.parse_args()
args_dict = args.__dict__

if not (args.ifd):
    sys.exit("ERROR: inputfastadir is required")

inputfastadir = os.path.abspath(args.ifd)

batch_out_path = args.o
fasta_files = [f for f in os.listdir(inputfastadir) if f.endswith('.fasta')]

if not len(fasta_files):
    sys.exit("ERROR: no FASTA files found")

run_resfinder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'run_resfinder.py')

arg_coll = []
for x, y in args_dict.items():
    if x not in ['ifd', 'o'] and y is not False and y is not None:
        if y is True:
            arg_coll.append(f'-{x}')
        else:
            arg_coll.append(f'-{x} {y}')

arg_str = ' '.join(arg_coll)

for f in fasta_files:
    inputfasta = os.path.join(inputfastadir, f)
    out_path = os.path.join(batch_out_path, f.replace('.fasta', ''))
    print(f'exec python3 {run_resfinder} -ifa {inputfasta} -o {out_path} {arg_str}')
    os.system(f'exec python3 {run_resfinder} -ifa {inputfasta} -o {out_path} {arg_str}')

time_taken = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
print(f'Processed {len(fasta_files)} FASTA files (Time taken: {time_taken})')

sys.exit()
