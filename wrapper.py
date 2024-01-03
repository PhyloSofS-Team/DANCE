#!/usr/bin/env python3
from src.python import mp_align, mp_cifAlignement, mp_cifConverter, mmseqs, mp_write_stats

def main():

    identity = 0.8  
    coverage = 0.8 
    bool_weight = True
    cif_dir = 'test/cifs/' # Directory containing CIF files
    mmseqs_dir = 'test/mmseqs_output/' # Directory to store the MMseqs2 output files
    aligned_dir = 'test/aligned/' # Directory to store the aligned sequences
    models_dir = 'test/models/' # Directory to store the models
    mf_name = 'test/output.fa' # Name of the multifasta file containing the all the sequences

    cif_alignment_options = {
        'c': True,  # Enable centermass and alignment on Calpha only
        'w': bool_weight,  # Enable weighted alignment
        'a': True,  # Enable output option aln
        'r': True,  # Enable output option RMSD
        'b': True,  # Enable output option raw coords
        'u': True,  # Enable output option for removed sequence information
        'p': True,  # Enable output option pdb
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
    mp_cifAlignement.run_cif_alignment(aligned_dir, cif_alignment_options)

    # Compute statistics of the models
    mp_write_stats.main(use_weights=bool_weight, directory=models_dir)

if __name__ == '__main__':
    main()