#!/usr/bin/env python3
from src.python import mp_align, mp_cifAlignement, mp_cifConverter, mmseqs, mp_write_stats

def main():
    options = {
        'c': True,  # Centermass et alignement sur Calpha uniquement
        'w': True,  # Alignement pondéré
        'a': True,  # Activer l'option de sortie aln
        'r': True,  # Activer l'option de sortie RMSD
        'b': True,  # Activer l'option de sortie raw coords
        'u': True,  # Activer l'option de sortie des informations de séquences retirées
        'p': True,  # Activer l'option de sortie pdb
        'd': 'test/cifs/',
        'o': 'test/models/'
    }

    mp_cifConverter.process_files('test/cifs/', 'test/output.fa')
    mmseqs.process_mmseqs('0.8', '0.8', 'test/output.fa', 'test/mmseqs_output')
    mp_align.process_files('test/output.fa', 'test/mmseqs_output/clusterDB_80_80.tsv', 'test/aligned')
    mp_cifAlignement.run_cif_alignment('test/aligned', options)
    mp_write_stats.main(use_weights=True, directory="test/models/")

if __name__ == '__main__':
    main()
