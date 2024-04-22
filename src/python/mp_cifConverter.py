import os
import subprocess
import multiprocessing
from tqdm import tqdm
import argparse


def setup_environment():
    env_path = os.path.join(os.environ["SCRIPT_DIR"], ".env")
    with open(env_path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value

    cif_converter_path_env = os.getenv("CIF_CONVERTER_PATH", "bin/cifConverter")

    if os.path.isabs(cif_converter_path_env):
        cif_converter_path = cif_converter_path_env
    else:
        cif_converter_path = os.path.join(
            os.environ["SCRIPT_DIR"], cif_converter_path_env
        )

    return cif_converter_path


def process_cif(file, cif_converter_args):
    command = [cif_converter_args[0]] + cif_converter_args[1:] + [file]
    try:
        result = subprocess.run(command, capture_output=True, text=True)
        return result.stdout
    except Exception as e:
        print(f"Error processing {file}: {e}")
        return ""


def process_cif_wrapper(args):
    return process_cif(*args)


def list_files(directory, extension):
    return [
        os.path.join(directory, file)
        for file in os.listdir(directory)
        if file.endswith(extension)
    ]


def process_files(
    input_path,
    output_file,
    num_workers=None,
    cif_converter_path=None,
    monomer_only=False,
):
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    if cif_converter_path is None:
        cif_converter_path = setup_environment()

    cif_converter_args = [cif_converter_path]
    if monomer_only:
        cif_converter_args.append("--monomer-only")

    print("Extracting sequences from CIF files:")
    cif_files = []
    if os.path.isfile(input_path) and input_path.endswith(".cif"):
        print(
            "Warning: Processing a single file. This script is optimized for parallel processing of multiple CIF files."
        )
        cif_files = [input_path]
    elif os.path.isdir(input_path):
        cif_files = list_files(input_path, ".cif")
    elif os.path.isfile(input_path):
        with open(input_path, "r") as file:
            cif_files = [line.strip() for line in file if line.strip()]

    file_args = [
        (file, cif_converter_args) for file in cif_files
    ]  # Pass the list of arguments
    with multiprocessing.Pool(processes=num_workers) as pool:
        sequences = list(
            tqdm(pool.imap(process_cif_wrapper, file_args), total=len(cif_files))
        )

    with open(output_file, "w") as f:
        for seq in sequences:
            f.write(seq)


def main():
    parser = argparse.ArgumentParser(description="Process CIF files.")

    parser.add_argument(
        "-i",
        "--input",
        help="A .cif file, a directory containing CIF files, or a file with a newline-separated list of CIF files",
    )
    parser.add_argument(
        "-m",
        "--cifConverterPath",
        required=True,
        help="Path to the cifConverter binary.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output.fa",
        help="Output file name (default: output.fa)",
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes (default: number of CPU cores)",
    )
    parser.add_argument(
        "--monomer-only", action="store_true", help="Only process monomeric assemblies."
    )
    args = parser.parse_args()
    process_files(
        args.input, args.output, args.workers, args.cifConverterPath, args.monomer_only
    )


if __name__ == "__main__":
    main()
