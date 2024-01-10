# DANCE: Dimensionality Analysis for protein Conformational Exploration

DANCE is designed to process a collection of protein 3D structures provided in Crystallographic Information File (CIF) format, producing protein or protein family-specific conformational ensembles in either CIF or PDB format. By leveraging amino acid sequence similarities, it effectively groups and aligns the inputted 3D structures to create these ensembles. Additionally, the tool compiles these ensembles and extracts principal components to define the 'linear motion manifold', a fundamental representation of conformational dynamics. It also quantifies the intrinsic dimensionality of these manifolds, offering insights into the structural variations and flexibility of proteins.

## Description of the Algorithm

DANCE's algorithm unfolds in the following steps:

- **Sequence Extraction:** Extracts amino acid sequences from CIF files, retaining the first model in cases of multiple models and including unresolved residues as lowercase letters or 'X' symbols. Sequences with less than 5 non-'X' residues are filtered out.

- **Sequence Clustering:** Utilizes MMseqs2 to cluster sequences based on user-defined similarity and coverage levels, defaulting to 80%. Outputs a TSV file detailing the clusters.

- **Multiple Sequence Alignments:** Aligns sequences within each cluster using MAFFT, removes columns with only 'X's or gaps, and orders sequences by their PDB codes for efficiency.

- **Structure Extraction and Generation of Conformational Ensembles:** This step involves extracting 3D coordinates of backbone atoms from CIF files, reconstructing missing atoms, and discarding incomplete residues or chains. It then groups conformations based on sequence clustering and superimposes them using the Quaternion Characteristic Polynomial method, reducing structural redundancy. The final output is a conformational ensemble in PDB/CIF format, accompanied by corresponding Multiple Sequence Alignments (MSAs) and Root Mean Square Deviation (RMSD) matrices.

- **Extraction of Linear Motions:** Performs Principal Component Analysis (PCA) on each ensemble's 3D coordinates to identify and interpret linear motions in the protein structures.


While each of these steps is managed by an individual Python script capable of functioning autonomously, the complete pipeline that invokes these scripts in sequence is implemented in wrapper.py. This approach allows users the flexibility to either utilize the entire pipeline for comprehensive analysis or employ individual scripts for specific tasks as needed.

### Getting Started

### Dependencies

- **Python 3.x**
- **Operating Systems:** Linux, macOS
- **External Tools:**
  - [Mafft](https://mafft.cbrc.jp/alignment/software/)
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
3. **Configure the Environment:**

- Open the `.env` file in a text editor.
- Modify the paths for `MAFFT_PATH` and `MMSEQS_PATH` to reflect the installation paths of Mafft and MMseqs2 on your system.

4. **Install Python Dependencies:**
```
pip install -r requirements.txt
```
5. **Compile the Binaries:**
```
cmake CMakeLists.txt
make
```
6. **Run the Exemple:**

You may want to run the pipeline on the small set of given cifs:
```
python wrapper.py
```