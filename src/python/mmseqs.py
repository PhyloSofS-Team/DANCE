import argparse
import subprocess
import os
import sys
from typing import List

env_path = '.env'

with open(env_path) as f:
    for line in f:
        if line.strip() and not line.startswith('#'):
            key, value = line.strip().split('=', 1)
            os.environ[key] = value

mmseqs_path = os.getenv('MMSEQS_PATH', 'mmseqs/bin/mmseqs')

def run_command(command: List[str]):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(e)

def rewrite_tsv(input_name: str, output_name: str):
    newfile = []
    newline = []
    
    with open(input_name, 'r') as f:
        for index, line in enumerate(f):
            elements = line.split()
            if elements[0] == elements[1]:
                if newline:  
                    newfile.append(sorted(newline))
                newline = [elements[1]]
            else:
                newline.append(elements[1])
    
    if newline:
        newfile.append(sorted(newline))
    
    newfile.sort()
    
    with open(output_name, 'w') as f:
        for line in newfile:
            f.write(' '.join(line) + '\n')


def process_mmseqs(id: str, cov: str, multifasta: str, output_dir: str, num_workers: int = None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    db_dir = os.path.join(output_dir, "db_files")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    db_file = os.path.join(db_dir, 'DB')
    covc = int(float(cov) * 100)
    idc = int(float(id) * 100)
    idcov = f"{idc}_{covc}"

    cluster_db = os.path.join(db_dir, f"clusterDB_{idcov}")
    tsv_file = os.path.join(db_dir, f"clusterDB_{idcov}.tsv")

    thread_option = f"--threads {num_workers}" if num_workers is not None else ""

    run_command(f"{mmseqs_path} createdb {multifasta} {db_file}")
    run_command(f"{mmseqs_path} cluster {db_file} {cluster_db} {db_dir}/tmp --cov-mode 0 -c {cov} --min-seq-id {id} {thread_option}")
    run_command(f"{mmseqs_path} createtsv {db_file} {db_file} {cluster_db} {tsv_file}")

    final_tsv = os.path.join(output_dir, f"clusterDB_{idcov}.tsv")
    rewrite_tsv(tsv_file, final_tsv)

def main():
    parser = argparse.ArgumentParser(description="Process FASTA files with mmseqs.")
    parser.add_argument('-i', '--id', required=True, help="Identity threshold")
    parser.add_argument('-c', '--cov', required=True, help="Coverage threshold")
    parser.add_argument('-m', '--multifasta', required=True, help="Path to the multifasta file")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for mmseqs files")
    parser.add_argument('-w', '--workers', type=int, help="Number of worker threads for mmseqs (optional)")
    args = parser.parse_args()

    process_mmseqs(args.id, args.cov, args.multifasta, args.output_dir, args.workers)

if __name__ == "__main__":
    main()