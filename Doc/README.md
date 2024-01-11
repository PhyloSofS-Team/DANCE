# Documentation

## CifConverter

CifConverter is a very simple tool written in C++ that extracts the amino acid sequence of a cif file in a multi-FASTA format.

#### Extraction of amino acid sequences:

- The code extracts amino acid sequences of all polypeptidic chains present in the input CIF files.
- It focuses on the **__atom_site.label_comp_id_** column to identify the names of residues with resolved 3D coordinates.

#### Inclusion of missing residues:

- Residues that are missing from the protein structure but are defined in the **__entity_poly_seq_** category are included in the sequence as lowercase letters.

#### Multiple models:

- In scenarios where multiple models are present in the CIF file, the code retains **only the first model** for sequence extraction.

#### Usage of the 'X' symbol:

- The 'X' symbol is utilized in the sequence in place of unknown amino acid types, modified amino acids without closely related natural counterparts, and residues missing from the structure and not defined in the __entity_poly_seq category_.

#### Sequence filtering:

- The code implements a filtering step where sequences comprising fewer than 5 non-'X' residues are excluded from further analysis.

### Standalone usage exemple:
#### In the terminal
```
$ bin/cifConverter test/cifs/1ake_final.cif 
>1AKEA	X-RAY DIFFRACTION	2.000000
MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
>1AKEB	X-RAY DIFFRACTION	2.000000
MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
```
#### To save the result in a fasta file
```
$ bin/cifConverter test/cifs/1ake_final.cif > 1ake_seq.fa
```

### Parallel launching using mp_cifConverter:

#### Whole directory:

To read the sequences of all the cifs in a directory and regroup them in one multi-FASTA file:
```
$ python src/python/mp_cifConverter.py test/cifs/ 
Extracting sequences from CIF files:
100%|████████████████████████████████████████| 792/792 [00:01<00:00, 686.70it/s]
```
The multi-FASTA file will be named  **output.fa** by default. You can use the **-o** option to chose your path and name.

#### Using a file list:

You can pass a file list as the main argument:
```
$ cat cifs.txt 
test/cifs/8afb_final.cif
test/cifs/8afc_final.cif
test/cifs/8afd_final.cif
test/cifs/8aht_final.cif
test/cifs/8aq5_final.cif
test/cifs/8aq7_final.cif
test/cifs/8azr_final.cif
test/cifs/8azv_final.cif
test/cifs/8azx_final.cif
test/cifs/8azy_final.cif
test/cifs/8b00_final.cif
test/cifs/8be6_final.cif
test/cifs/8be7_final.cif
test/cifs/8be8_final.cif
test/cifs/8be9_final.cif
test/cifs/8bea_final.cif
test/cifs/8bqf_final.cif
test/cifs/8cx5_final.cif
test/cifs/8dni_final.cif
test/cifs/8dnj_final.cif
test/cifs/8dnk_final.cif
test/cifs/8ebz_final.cif
test/cifs/8epw_final.cif
test/cifs/8h7f_final.cif
test/cifs/8h7h_final.cif
test/cifs/8onv_final.cif

$ python src/python/mp_cifConverter.py cifs.txt
Extracting sequences from CIF files:
100%|██████████████████████████████████████████████████| 26/26 [00:00<00:00, 496.12it/s]
```

## Mmseqs launcher

The goal of this script is to launch the [MMseqs2](https://github.com/soedinglab/MMseqs2) clustering on the extracted sequence with the desired identity threshold and coverage.

The coverage is bidirectional by default.

It also rewrites TSV file containing the ids of the created clusters to a one cluster per line file. 

#### Exemple usage:

```
python src/python/mmseqs.py -i 0.8 -c 0.8 -m test/output.fa -o test/mmseq_dir/
```
Where:
- **-i** is the identity threshold.
- **-c** is the coverage threshold.
- **-o** is the output directory for the Mmseqs2 files. The TSV will be outputted in this directory.

## Multiple sequence alignment:

This script first creates unaligned multi-FASTA file then perform multiple sequence alignment within each cluster using [MAFFT](https://mafft.cbrc.jp/alignment/software/), launched in parallel.

#### Exemple usage:

```
python src/python/mp_align.py -m test/output.fa -t test/mmseqs_output/clusterDB_80_80.tsv -d test/aligned_seq/
```
Where :
- **-m** is the path of the multifasta containing all the sequences.
- **-t** is the path of the TSV file containing the clusters.
- **-d** is the output directory for the aligned multi-FASTA files.

## CifAlignment

