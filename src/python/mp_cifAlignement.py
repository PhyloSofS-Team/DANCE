import sys
import glob
import subprocess
import time
import multiprocessing
import tqdm
import argparse
import logging
import os
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
    cmd = [cif_alignment_path, '-i', mf]
    for key, value in options.items():
        if value:
            if key in ['s', 'n', 'x', 'y', 'd', 'o']: 
                cmd.extend(['-' + key, str(value)])
            else:
                cmd.append('-' + key)
            
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    log_entry = f"Executing: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}\n"
    
    # Logging the command execution details
    logging.info(f"Executing: {' '.join(cmd)}")
    logging.info(f"STDOUT: {result.stdout}")
    if result.stderr:
        logging.error(f"STDERR: {result.stderr}")
    
    return log_entry

def process_cif_alignment(mf_files, options, log_filename=None):

    output_dir = options.get('o', '.')
    os.makedirs(output_dir, exist_ok=True)   
    if log_filename is None:
        current_time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        log_filename = f"cifAlignment-{current_time}.log"

    log_path = os.path.join(output_dir, log_filename)  
    logging.basicConfig(filename=log_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('Started logging.')

    logging.info('Launching cifAlignment')
    print('Creating aligned multistructure files:')
    t1 = time.time()
    with multiprocessing.Pool() as pool:
        log_entries = list(tqdm.tqdm(pool.imap_unordered(launchCifAlignmentWrapper, [(mf, options) for mf in mf_files]), total=len(mf_files)))
    t2 = time.time()

    logging.info(f"cifAlignment time: {t2 - t1}")
    logging.info('Finished logging.')

    return log_entries


def run_cif_alignment(mfdir, options, listfile=None):
    if listfile:
        with open(listfile, 'r') as f:
            mf_files = [line.strip() for line in f.readlines()]
    else:
        # Générer une liste de fichiers .fa, en excluant les fichiers *_unaligned.fa
        all_fa_files = set(glob.glob(os.path.join(mfdir, '*.fa')))
        unaligned_fa_files = set(glob.glob(os.path.join(mfdir, '*_unaligned.fa')))
        mf_files = list(all_fa_files - unaligned_fa_files)

    process_cif_alignment(mf_files, options)

def main():
    parser = argparse.ArgumentParser(description="Run cifAlignment with specified options.")
    parser.add_argument("-i", "--inputAln", help="Path to the aln file.")
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
    if args.inputAln:
        mf_files = [args.inputAln]
    elif args.listfile:
        with open(args.listfile, 'r') as f:
            mf_files = [line.strip() for line in f.readlines()]
    else:
        raise ValueError("Either --inputAln or --listfile must be provided")

    run_cif_alignment(mf_files, options)


if __name__ == "__main__":
    main()