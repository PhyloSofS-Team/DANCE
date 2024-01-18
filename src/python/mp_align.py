import sys
import tqdm
import subprocess
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
            fasta_data = ('>' + multi_fasta_file.read()).split('\n>')
            
        
        Dict = {}
        for data in fasta_data:
            if data:
                header, sequence = data.split('\n', 1)
                seq_id = header.split()[0].strip('>') 
                Dict[seq_id] = sequence.strip()

        lclusters_unaligned = []

        print('Reading members and creating unaligned multifasta files:')
        for cluster in tqdm.tqdm(clusters):
            if cluster:
                cluster_id = cluster[0]
                namefile_unaligned = os.path.join(name_dr, f'{cluster_id}_unaligned.fa')
                with open(namefile_unaligned, 'w') as file:
                    for seq_id in cluster:
                        sequence = Dict.get(seq_id, "")
                        if sequence:
                            file.write(f'>{seq_id}\n{sequence}\n')
                    if len(cluster) > 1:
                        lclusters_unaligned.append(namefile_unaligned)
        
        return lclusters_unaligned

    except Exception as e:
        print(f"Error occurred: {e}")

def launch_mafft_on_clusters(lclusters_unaligned, name_dr, num_workers=None):
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()  

    os.makedirs(name_dr, exist_ok=True)
    lclusters_aligned = [cluster.replace('_unaligned', '') for cluster in lclusters_unaligned]

    print('Launching Mafft on multifasta files:')

    try:
        with multiprocessing.Pool(processes=num_workers) as pool:
            list(tqdm.tqdm(pool.imap_unordered(launchMafft, zip(lclusters_unaligned, lclusters_aligned)), total=len(lclusters_aligned)))
    except Exception as e:
        sys.stderr.write(f"Error in launch_mafft_on_clusters: {e}\n")


def process_files(multifasta, tsv, name_dr, num_workers=None):
    lclusters_unaligned = create_unaligned_fasta_files(multifasta, tsv, name_dr)
    launch_mafft_on_clusters(lclusters_unaligned, name_dr, num_workers)


def main():
    parser = argparse.ArgumentParser(description="Run MAFFT on multiple FASTA files.")
    parser.add_argument('-m', '--multifasta', required=True, help="Path to the multifasta file")
    parser.add_argument('-t', '--tsv', required=True, help="Path to the TSV file containing clusters information")
    parser.add_argument('-d', '--output_dir', required=True, help="Output directory for aligned FASTA files")
    parser.add_argument('-w', '--workers', type=int, default=None, help="Number of worker processes (default: number of CPU cores)")

    args = parser.parse_args()

    process_files(args.multifasta, args.tsv, args.output_dir, args.workers)

if __name__ == "__main__":
    main()