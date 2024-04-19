import numpy as np
from Bio.SeqUtils import seq1


def read_pdb_info(pdb_file):
    """
    Read the header of a pdb file and extract sequence length, number of conformations,
    and names of the conformations.
    Arguments:
    pdb_file -- The path to the pdb file.
    Returns:
    A tuple containing sequence length, number of conformations, and names of the conformations.
    """
    with open(pdb_file) as f:
        line_1 = f.readline()
        if line_1.split()[0] != "NUMMDL":
            return (None, None, None)
        else:
            num_mdl = int(line_1.split()[1])
            line_2 = f.readline()
            seq_len = int(line_2.split()[1])
            chains_names = [f.readline().split()[2] for _ in range(num_mdl)]
            return (seq_len, num_mdl, chains_names)


def parse_pdb_file(
    pdb_file_path,
    CA_only=False,
    output_gaps=False,
    output_alignment=False,
    trim_single_residue_columns=False,
    output_chain_names=False,
    output_sequences=False,
):
    """
    Parses a _mm.pdb file and extracts relevant information.

    Args:
        pdb_file_path (str): The path to the _mm.pdb file.
        CA_only (bool): If True, extracts only CA (carbon alpha) atoms.
        output_gaps (bool): If True, outputs the gaps matrix.
        output_alignment (bool): If True, outputs the alignment matrix.
        trim_single_residue_columns (bool): If True, removes columns containing only one residue.
        output_chain_names (bool): If True, outputs the names of the chains.
        output_sequences (bool): If True, outputs the sequences.
    Returns:
        dict: A dictionary containing the parsed information.
    """

    # Read basic PDB file info
    sequence_length, num_models, chain_names = read_pdb_info(pdb_file_path)
    if sequence_length is None:
        raise ValueError("PDB multi-model information is missing!")

    floats_per_residue, atoms_per_residue = (3, 1) if CA_only else (12, 4)

    parsed_info = {}
    coordinates = np.zeros((sequence_length * floats_per_residue, num_models))
    gaps = np.zeros((sequence_length * floats_per_residue, num_models))

    if output_alignment:
        alignment = np.full((num_models, sequence_length), "-")

    model_index = 0
    residue_id_0 = -1

    with open(pdb_file_path) as pdb_file:
        for line in pdb_file:
            if line.startswith(("ATOM", "HETATM")):

                parsed_line = [
                    line[i:j]
                    for i, j in [
                        (0, 6),
                        (6, 11),
                        (12, 16),
                        (17, 20),
                        (21, 22),
                        (22, 26),
                        (30, 38),
                        (38, 46),
                        (46, 54),
                    ]
                ]
                if CA_only and parsed_line[2].strip() != "CA":
                    continue

                residue_id = int(parsed_line[5]) - 1
                if output_alignment and (residue_id != residue_id_0):
                    alignment[model_index, residue_id] = seq1(parsed_line[3]).upper()
                    residue_id_0 = residue_id

                atom_id = int(parsed_line[1]) - 1
                coord_index = (residue_id * floats_per_residue) + (
                    (atom_id % atoms_per_residue) * 3
                )

                coordinates[coord_index : coord_index + 3, model_index] = np.array(
                    parsed_line[6:9], dtype=float
                )
                gaps[coord_index : coord_index + 3, model_index] = 1

            elif line.startswith("MODEL"):
                model_index = int(line.split()[1]) - 1

    if trim_single_residue_columns:
        columns_to_keep = np.where(gaps.sum(axis=1) > 1)[0]
        coordinates = coordinates[columns_to_keep]
        gaps = gaps[columns_to_keep]
        if output_alignment:
            alignment = alignment[:, list(set(columns_to_keep // floats_per_residue))]

    if output_chain_names:
        parsed_info["chain_names"] = chain_names
    if output_alignment:
        parsed_info["alignment"] = alignment.T
    if output_gaps:
        parsed_info["gaps"] = gaps
    if output_sequences:
        if output_alignment is False:
            print(
                "Warning: output_alignment is False, sequences will not be outputted!"
            )
        else:
            sequences = ["".join(seq).replace("-", "") for seq in alignment]
            parsed_info["sequences"] = sequences

    parsed_info["coordinates"] = coordinates

    return parsed_info


# Example:
# dic = parse_pdb_file('test/models/1AKEA_1AKEA_mm.pdb', output_alignment=True, CA_only=False,output_sequences=True)
