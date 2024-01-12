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

#### Using a file list:

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

## cifAlignment

### Overview

cifAlignment superimposes structures using a multi-FASTA file. It first chooses a reference structure, then aligns other structures to it by computing the optimal least-squares rotation matrix. This process minimizes the Root Mean Square Deviation (RMSD) between each conformation and the reference. The output is provided in PDB or CIF file formats. Additionally, cifAlignment can generate a corresponding Multiple Sequence Alignment (MSA) in multi-FASTA format and an RMSD matrix of the structures.

#### Choice of the reference

We select the reference conformation for superimposition based on the amino acid sequence that is most representative of the MSA. To do this, we first determine the consensus sequence $s^{\*}$ by identifying the most frequent symbol at each position. "X" symbols are considered equivalent to gaps. Therefore, each position is described by a 21-dimensional vector representing the frequencies of the 20 amino acid types and the gaps. In case of ambiguity, we prefer an amino acid over a gap, and a more frequent amino acid over a less frequent one. Then, we compute a score for each sequence $s$ in the MSA reflecting its similarity to $s^{\*}$, expressed as:

$score(s) = ∑_{i=1}^P\sigma(s_i, s^{\*}_i)$,

where $P$ is the number of positions in the MSA and $\sigma(s_i, s^{\*}_i)$  is the substitution score between the amino acid $s_i$ at position i in sequence s and the consensus symbol $s^{\*}_i$ at the same position. We use the substitution matrix BLOSUM62, setting the gap score to $min_{a,b}(\sigma(a,b))-1=-5$.

The choice of reference can be forced by using the **--numReferences** parameter.

#### Removal of isolated residues

cifAlignment has the ability to remove residues that are too isolated in a sequence. The idea behind this function is that an isolated residue in a sequence is more likely to be misaligned. This feature uses two parameters, **--continentSize** and **--isolationDistance**.

**--continentSize** defines the size from which a contiguous section of residue, that is, a part of the sequence containing neither gap nor "X", is considered as a continent and can never be removed from a sequence. Residue sections smaller than **--continentSize** are considered as islands.

**--isolationDistance** defines the distance from which an island is removed from the sequence. This occurs when the island in question is located more than **--isolationDistance** positions from the nearest continent.

If one wants to ensure that no residues are removed from the sequence, then either of the two parameters **--continentSize** or **--isolationDistance** must be set to 0.

#### Removal of empty columns 

If a column contains only gaps or 'X's, it is removed from the alignment.

#### Structural alignment on the reference

The structural alignment of the conformations to the reference is done by minimizing the Root Mean Squared Deviation (RMSD) between the common set of occupied positions in the sequence alignment between the reference and the considered conformation

##### Calpha or whole backbone

The alignment can be done by using only the alpha carbon with the option **--calpha**, or with all the backbone atoms without this option.

##### Removal of redundant structures

This step also computes the values of the RMSD matrix, allowing removing redundant structures by using a RMSD threshold defined by **--similarity**.
If two structures have less than **--similarity** RMSD, and identical amino acid sequence, we **remove the last in the alphabetical order**.
If two structures have less than **--similarity** RMSD but one is subset of the other, we **remove the shortest conformation**.
If two structures have less than **--similarity** RMSD but have different sets of residues, we **keep them both**.

##### Minimal substet in common for structural alignment

If two structures have less than **--commonResAln** residues in common, their RMSD is set to NAN. If a structure has **NAN RMSD with the reference**, it is **removed from the ensemble**.

##### Weighted structural alignment

cifAlignment offers the possibility to weight the structural alignment by giving relative wheights to the common subset of residues between the two aligned structures. We derive those weights from the global coverage of the position in the MSA. To activate this option, use the **--weighted** option.

### Output

cifAlignment saves the resulting conformational ensemble as a multi-model file in either PDB or CIF format. It's important to note that these models can display different amino acid sequences. Additionally, cifAlignment generates and outputs the corresponding multiple sequence alignments (MSA) in FASTA format, alongside a matrix of all-to-all pairwise RMSDs.

#### PDB file

The output is a PDB file containing multiple models. The header of the file contains information about the number of models and the length of the alignment. It also contains a list of the ids of the different models of the PDB. Please note that this PDB file is non-standard because the amino acids present in the models may differ.
We use the residue sequence number at columns 23 to 26 of the PDB file in order to represent the position in the alignment of each residue.
```
$ head test/models/1AKEA_1AKEA_mm.pdb 
NUMMDL    35                                                                    
SEQLEN    215
1 : 1AKEA
2 : 1AKEB
3 : 1E4YA
4 : 1E4YB
5 : 3HPQA
6 : 3HPQB
7 : 3HPRA
8 : 3HPRB
...
```