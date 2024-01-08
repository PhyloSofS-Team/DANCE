import sys
import glob
import subprocess
import time
import multiprocessing
import tqdm
import argparse
import os
import logging
import datetime

env_path = '.env'

with open(env_path) as f:
    for line in f:
        if line.strip() and not line.startswith('#'):
            key, value = line.strip().split('=', 1)
            os.environ[key] = value

cif_alignment_path = os.getenv('CIF_ALIGNMENT_PATH', 'bin/cifAlignment')

def launchCifAlignmentWrapper(args):
    return launchCifAlignment(*args)

def launchCifAlignment(mf, options):
    mf_base = os.path.splitext(os.path.basename(mf))[0]
    log_filename = f"log_{mf_base}.log"
    output_dir = options.get('o', '.') 
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, log_filename)  
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_path)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    cmd = [cif_alignment_path, '-i', mf]
    for key, value in options.items():
        if value:
            if key in ['s', 'n', 'x', 'y', 'd', 'o']:
                cmd.extend(['-' + key, str(value)])
            else:
                cmd.append('-' + key)

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    logger.info(f"Executing: {' '.join(cmd)}")
    logger.info(f"STDOUT: {result.stdout}")
    if result.stderr:
        logger.error(f"STDERR: {result.stderr}")

    logger.removeHandler(handler)

    return result.stdout, result.stderr



def process_cif_alignment(mf_files, options):
    print('Creating aligned multistructure files:')
    t1 = time.time()
    with multiprocessing.Pool() as pool:
        for _ in tqdm.tqdm(pool.imap_unordered(launchCifAlignmentWrapper, [(mf, options) for mf in mf_files]), total=len(mf_files)):
            pass
    t2 = time.time()

    


def run_cif_alignment(mfdir, options, listfile=None):
    if listfile:
        with open(listfile, 'r') as f:
            mf_files = [line.strip() for line in f.readlines()]
    else:
        all_fa_files = set(glob.glob(os.path.join(mfdir, '*.fa')))
        unaligned_fa_files = set(glob.glob(os.path.join(mfdir, '*_unaligned.fa')))
        mf_files = list(all_fa_files - unaligned_fa_files)

    process_cif_alignment(mf_files, options)

def main():
    parser = argparse.ArgumentParser(description="Run cifAlignment with specified options.")
    parser.add_argument("-i", "--alnDir", help="Path to the aln directory.")
    parser.add_argument("-d", "--cifDir", help="Path to the cif directory.")
    parser.add_argument("-o", "--outputDir", help="Path to the output directory.")
    parser.add_argument("-c", "--centermass", action="store_true", help="Centermass and alignment on Calpha only.")
    parser.add_argument("-w", "--weighted", action="store_true", help="Enable weighted alignment.")
    parser.add_argument("-s", "--similarity", type=float, help="Set RMSD similarity threshold for conformation removing.")
    parser.add_argument("-p", "--outputPdb", action="store_true", help="Enable output pdb file.")
    parser.add_argument("-f", "--outputCif", action="store_true", help="Enable output cif file.")
    parser.add_argument("-a", "--outputAln", action="store_true", help="Enable output aln file.")
    parser.add_argument("-r", "--outputRmsd", action="store_true", help="Enable output RMSD file option.")
    parser.add_argument("-u", "--outputRemoved", action="store_true", help="Enable output of removed sequences file information.")
    parser.add_argument("-b", "--outputRawCoords", action="store_true", help="Enable output raw coords file option.")
    parser.add_argument("-n", "--numReferences", type=int, help="Set the number of references.")
    parser.add_argument("-x", "--continentSize", type=int, help="Set the continent size (strictly superior to).")
    parser.add_argument("-y", "--isolationDistance", type=int, help="Set the isolation distance (superior or equal to).")
    parser.add_argument("-l", "--listfile", help="Path to the file containing a list of .fa files.")
    args = parser.parse_args()

    options = {
        'c': args.centermass,
        'w': args.weighted,
        'p': args.outputPdb,
        'f': args.outputCif,
        'a': args.outputAln,
        'r': args.outputRmsd,
        'u': args.outputRemoved,
        'b': args.outputRawCoords,
        's': args.similarity,
        'n': args.numReferences,
        'x': args.continentSize, 
        'y': args.isolationDistance,
        'd': args.cifDir,  
        'o': args.outputDir  
    }
    
    if not args.alnDir and not args.listfile:
        raise ValueError("Either --alnDir or --listfile must be provided, use -h for help.")
    else:
        run_cif_alignment(args.alnDir, options, args.listfile)


if __name__ == "__main__":
    main()