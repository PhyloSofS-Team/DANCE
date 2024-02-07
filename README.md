# DANCE: Dimensionality Analysis for protein Conformational Exploration

DANCE is designed to process a set of input protein 3D structures provided in Crystallographic Information File (CIF) format and output protein or protein family-specific conformational collections in CIF or PDB format. By leveraging amino acid sequence similarities, it effectively groups and aligns the inputted 3D structures to create these collections (or ensembles). Additionally, DANCE extracts the principal components from each conformational collection, using classical Principal Component Analysis. The principal components represent directions in the 3D space for every atom and can be interpreted as 'linear motions'. The eigenvalue associated to each component reflects its contribution to the total positional variance of the collection.

## Description of the Algorithm

DANCE's algorithm unfolds in the following steps:

- **Sequence Extraction:** Extracts amino acid sequences from CIF files, retaining the first model in cases of multiple models and including unresolved residues as lowercase letters or 'X' symbols. Sequences with less than 5 non-'X' residues are filtered out.

- **Sequence Clustering:** Utilizes MMseqs2 to cluster sequences based on user-defined similarity and coverage levels, defaulting to 80%. Outputs a TSV file detailing the clusters.

- **Multiple Sequence Alignments:** Aligns sequences within each cluster using MAFFT, and orders sequences by their PDB codes for efficiency.

- **Structure Extraction and Generation of Conformational Ensembles:** This step involves extracting 3D coordinates of backbone atoms from CIF files, reconstructing missing atoms, and discarding incomplete residues or chains. It then groups conformations based on sequence clustering and superimposes them using the Quaternion Characteristic Polynomial method, reducing structural redundancy. The final output is a conformational ensemble in PDB/CIF format, accompanied by corresponding Multiple Sequence Alignments (MSAs) and Root Mean Square Deviation (RMSD) matrices.

- **Extraction of Linear Motions:** Performs Principal Component Analysis (PCA) on each ensemble's 3D coordinates to identify and interpret linear motions in the protein structures. The output is a CSV file containing summary statistics for the created ensembles.


While each of these steps is managed by an individual Python script capable of functioning autonomously, the complete pipeline that invokes these scripts in sequence is implemented in wrapper.py. This approach allows users the flexibility to either utilize the entire pipeline for comprehensive analysis or employ individual scripts for specific tasks as needed.

### Getting Started

### Dependencies

- **Python 3.x**
- **Operating Systems:** Linux, macOS
- **External Tools:**
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)
  - [MMseqs2](https://github.com/soedinglab/MMseqs2)
- **Python Libraries:** Installed via `requirements.txt`.

### Installing

1. **Clone the Repository:**
```
git clone https://github.com/PhyloSofS-Team/DANCE
```

2. **Navigate to the DANCE Directory:**
```
cd DANCE
```  
3. **Install Python Dependencies:**
```
pip install -r requirements.txt
```  
4. **Compile the Binaries:**
```
mkdir build
cd build
cmake ../
make
```
5. **Configure the Environment:**

- Open the `.env` file in a text editor and set the path variables `MAFFT_PATH` and `MMSEQS_PATH` to the installation paths of Mafft and MMseqs2 on your system.

6. **Run the Exemple:**

You may want to run the pipeline on the small set of given cifs:
```
$ python wrapper.py
Extracting sequences from CIF files:
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 792/792 [00:01<00:00, 632.18it/s]
[...]
Reading members and creating unaligned multifasta files:
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 140/140 [00:00<00:00, 22793.36it/s]
Launching Mafft on multifasta files:
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 78/78 [00:16<00:00,  4.71it/s]
Creating aligned multistructure files:
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 78/78 [00:06<00:00, 12.94it/s]
Computing statistics:
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████| 78/78 [00:01<00:00, 57.74it/s]
```
