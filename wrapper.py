#!/usr/bin/env python3
from src.python import mp_align, mp_cifAlignment, mp_cifConverter, mmseqs, mp_write_stats

def main():

    identity = 0.8  
    coverage = 0.8 
    bool_weight = False # Enable weighted structural alignment. Each position is weighted by it's coverage in the alignment.
    cif_dir = 'test/cifs/' # Directory containing CIF files
    mmseqs_dir = 'test/mmseqs_output/' # Directory to store the MMseqs2 output files
    aligned_dir = 'test/aligned/' # Directory to store the aligned multifasta files of the ensembles before model building
    models_dir = 'test/models/' # Directory to store the 3D models 
    mf_name = 'test/output.fa' # Name of the multifasta file containing all the sequences extracted from the CIF files of cif_dir

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
    mp_cifConverter.process_files(cif_dir, mf_name)

    # Create ensembles of sequences with MMseqs2
    mmseqs.process_mmseqs(identity, coverage, mf_name, mmseqs_dir)
  
    # Align sequences of the ensembles with Mafft
    mp_align.process_files(mf_name, f'{mmseqs_dir}{cluster_db_filename}', aligned_dir)

    # Run CIF alignment with specified options to create the 3d models
    mp_cifAlignment.run_cif_alignment(aligned_dir, cif_alignment_options)

    # Compute statistics of the models
    mp_write_stats.main(use_weights=bool_weight, directory=models_dir) # to write the stats, the options a, r and b must be enabled in cif_alignment_options

if __name__ == '__main__':
    main()