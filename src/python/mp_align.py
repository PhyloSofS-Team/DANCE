import sys
import tqdm
import subprocess
import multiprocessing
import argparse
import os


def setup_environment():
    env_path = os.path.join(os.environ["SCRIPT_DIR"], ".env")
    with open(env_path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value

    mafft_path_env = os.getenv("MAFFT_PATH", "mafft/bin/mafft")

    if os.path.isabs(mafft_path_env):
        mafft_path = mafft_path_env
    else:
        mafft_path = os.path.join(os.environ["SCRIPT_DIR"], mafft_path_env)

    return mafft_path


def launchMafft(args):
    try:
        cluster_unaligned, cluster_aligned, mafft_path = args
        with open(cluster_aligned, "w") as outfile:
            subprocess.run(
                [
                    mafft_path,
                    "--quiet",
                    "--auto",
                    "--amino",
                    "--preservecase",
                    cluster_unaligned,
                ],
                stdout=outfile,
            )
    except Exception as e:
        sys.stderr.write(f"Error in launchMafft: {e}\n")


def create_unaligned_fasta_files(mf, tsv, name_dr):
    os.makedirs(name_dr, exist_ok=True)
    try:
        with open(tsv, "r") as tsv_file:
            clusters = [line.strip().split() for line in tsv_file if line.strip()]

        with open(mf, "r") as multi_fasta_file:
            fasta_data = (">" + multi_fasta_file.read()).split("\n>")

        Dict = {}
        for data in fasta_data:
            if data:
                header, sequence = data.split("\n", 1)
                seq_id = header.split()[0].strip(">")
                Dict[seq_id] = sequence.strip()

        lclusters_unaligned = []

        print("Reading members and creating unaligned multifasta files:")
        for cluster in tqdm.tqdm(clusters):
            if cluster:
                cluster_id = cluster[0]
                namefile_unaligned = os.path.join(name_dr, f"{cluster_id}_unaligned.fa")
                with open(namefile_unaligned, "w") as file:
                    for seq_id in cluster:
                        sequence = Dict.get(seq_id, "")
                        if sequence:
                            file.write(f">{seq_id}\n{sequence}\n")
                    if len(cluster) > 1:
                        lclusters_unaligned.append(namefile_unaligned)

        return lclusters_unaligned

    except Exception as e:
        print(f"Error occurred: {e}")


def launch_mafft_on_clusters(
    lclusters_unaligned, name_dr, num_workers=None, mafft_path=None
):

    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    if mafft_path is None:
        mafft_path = setup_environment()

    os.makedirs(name_dr, exist_ok=True)
    lclusters_aligned = [
        cluster.replace("_unaligned", "") for cluster in lclusters_unaligned
    ]

    print("Launching Mafft on multifasta files:")

    mafft_args = [
        (cluster_unaligned, cluster_aligned, mafft_path)
        for cluster_unaligned, cluster_aligned in zip(
            lclusters_unaligned, lclusters_aligned
        )
    ]

    try:
        with multiprocessing.Pool(processes=num_workers) as pool:
            list(
                tqdm.tqdm(
                    pool.imap_unordered(launchMafft, mafft_args),
                    total=len(lclusters_aligned),
                )
            )
    except Exception as e:
        sys.stderr.write(f"Error in launch_mafft_on_clusters: {e}\n")


def process_files(multifasta, tsv, name_dr, num_workers=None, mafft_path=None):
    lclusters_unaligned = create_unaligned_fasta_files(multifasta, tsv, name_dr)
    launch_mafft_on_clusters(lclusters_unaligned, name_dr, num_workers, mafft_path)


def main():
    parser = argparse.ArgumentParser(description="Run MAFFT on multiple FASTA files.")
    parser.add_argument(
        "-m", "--mafft_path", required=True, help="Path to the MAFFT binary"
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Path to the multifasta file"
    )
    parser.add_argument(
        "-t",
        "--tsv",
        required=True,
        help="Path to the TSV file containing clusters information",
    )
    parser.add_argument(
        "-d",
        "--output_dir",
        required=True,
        help="Output directory for aligned FASTA files",
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes (default: number of CPU cores)",
    )

    args = parser.parse_args()

    process_files(args.fasta, args.tsv, args.output_dir, args.workers, args.mafft_path)


if __name__ == "__main__":
    main()
