import sys
import tqdm
import subprocess
import time
import multiprocessing
import argparse
import os

env_path = '.env'

with open(env_path) as f:
    for line in f:
        if line.strip() and not line.startswith('#'):
            key, value = line.strip().split('=', 1)
            os.environ[key] = value

mafft_path = os.getenv('MAFFT_PATH', 'mafft/bin/mafft')

def launchMafft(zip):
    try:
        with open(zip[1],'w') as outfile:
            subprocess.run([mafft_path, "--quiet", "--auto", "--amino", "--preservecase", zip[0]], stdout=outfile)
    except Exception as e:
        sys.stderr.write(f"Error in launchMafft: {e}\n")

def create_unaligned_fasta_files(mf, tsv, name_dr):
    os.makedirs(name_dr, exist_ok=True)
    try:
        with open(tsv, 'r') as tsv_file:
            clusters = [line.strip().split() for line in tsv_file if line.strip()]

        with open(mf, 'r') as multi_fasta_file:
            fasta_data = multi_fasta_file.read().split('\n>')
            
        Dict = {data.split('\n', 1)[0].split(sep='\t')[0].strip(): data.split('\n', 1)[1] for data in fasta_data}
        lclusters_unaligned = []

        print('Reading members and creating unaligned multifasta files:')
        for cluster in tqdm.tqdm(clusters):
            if cluster:
                cluster_id = cluster[0]
                namefile_unaligned = os.path.join(name_dr, f'{cluster_id}_unaligned.fa')
                with open(namefile_unaligned, 'w') as file:
                    for seq_id in cluster:
                        file.write(f'>{seq_id}\n{Dict.get(seq_id, "")}\n')
                    if len(cluster) > 1:
                        lclusters_unaligned.append(namefile_unaligned)
        
        return lclusters_unaligned

    except FileNotFoundError as e:
        sys.stderr.write(f"File not found: {e}\n")
    except Exception as e:
        sys.stderr.write(f"Error in create_unaligned_fasta_files: {e}\n")

def launch_mafft_on_clusters(lclusters_unaligned, name_dr):
    os.makedirs(name_dr, exist_ok=True)
    lclusters_aligned = [cluster.replace('_unaligned', '') for cluster in lclusters_unaligned]

    print('Launching Mafft on multifasta files:')
    t1 = time.time()
    try:
        with multiprocessing.Pool() as pool:
            list(tqdm.tqdm(pool.imap_unordered(launchMafft, zip(lclusters_unaligned, lclusters_aligned)), total=len(lclusters_aligned)))
    except Exception as e:
        sys.stderr.write(f"Error in launch_mafft_on_clusters: {e}\n")
    t2 = time.time()
    #print(f"{name_dr} alignment time: {t2 - t1}")

def process_files(mf, tsv, name_dr):
    lclusters_unaligned = create_unaligned_fasta_files(mf, tsv, name_dr)
    launch_mafft_on_clusters(lclusters_unaligned, name_dr)

def main():
    parser = argparse.ArgumentParser(description="Run MAFFT on multiple FASTA files.")
    parser.add_argument('-m', '--multifasta', required=True, help="Path to the multifasta file")
    parser.add_argument('-t', '--tsv', required=True, help="Path to the TSV file containing clusters information")
    parser.add_argument('-d', '--output_dir', required=True, help="Output directory for aligned FASTA files")
    
    args = parser.parse_args()

    process_files(args.multifasta, args.tsv, args.output_dir)

if __name__ == "__main__":
    main()
