#!/usr/bin/env python3script_directory = os.path.dirname(os.path.abspath(__file__))
import os
script_directory = os.path.dirname(os.path.abspath(__file__))
os.environ['SCRIPT_DIR'] = script_directory

from src.python import mp_align, mp_cifAlignment, mp_cifConverter, mmseqs, mp_write_stats

def main():

    identity = 0.8  
    coverage = 0.8 
    bool_weight = False # Enable weighted structural alignment. Each position is weighted by it's coverage in the alignment.

    # You can change the following paths to your own directories:
    cif_dir = os.path.join(script_directory,'test/cifs/') # Directory containing CIF files
    mmseqs_dir = os.path.join(script_directory,'test/mmseqs_output/') # Directory to store the MMseqs2 output files
    aligned_dir = os.path.join(script_directory,'test/aligned/') # Directory to store the aligned multifasta files of the ensembles before model building
    models_dir = os.path.join(script_directory,'test/models/') # Directory to store the 3D models 
    mf_name = os.path.join(script_directory,'test/output.fa') # Name of the multifasta file containing all the sequences extracted from the CIF files of cif_dir
    num_workers = None # Number of worker processes (default: number of CPU cores)
    cif_alignment_options = {
        'c': False,  # Alignment and center of mass calculation is done only with the CA atoms
        'w': bool_weight,  # Enable weighted alignment
        'a': True,  # Enable the output of the multifasta file of the aligned sequences corresponding to the models
        'r': True,  # Enable output option RMSD matrix of the conformations corresponding to the models
        'b': True,  # Enable output option raw coordinates in binary format of the conformations corresponding to the models
        'u': True,  # Enable output option for removed sequence information
        'p': True,  # Enable output of the model in a PDB file
        'f': False,  # Enable output of the model in a CIF file
        'n': 1,  # Number of models to build
        'd': cif_dir,
        'o': models_dir
    }

    cluster_db_file_suffix = f"{int(identity * 100)}_{int(coverage * 100)}"
    cluster_db_filename = f"clusterDB_{cluster_db_file_suffix}.tsv"

    # Extract sequences from CIF files
    mp_cifConverter.process_files(input_path=cif_dir, output_file=mf_name, num_workers=num_workers)

    # Create ensembles of sequences with MMseqs2
    mmseqs.process_mmseqs(id=identity, cov=coverage, multifasta=mf_name, output_dir=mmseqs_dir, num_workers=num_workers)
  
    # Align sequences of the ensembles with Mafft
    mp_align.process_files(multifasta=mf_name, tsv=f'{mmseqs_dir}{cluster_db_filename}', name_dr=aligned_dir, num_workers=num_workers)

    # Run CIF alignment with specified options to create the 3d models
    mp_cifAlignment.run_cif_alignment(mfdir=aligned_dir, options=cif_alignment_options, num_workers=num_workers)

    # Compute statistics of the models
    mp_write_stats.main(use_weights=bool_weight, directory=models_dir, num_workers=num_workers) # to write the stats, the options a, r and b must be enabled in cif_alignment_options

if __name__ == '__main__':
    main()