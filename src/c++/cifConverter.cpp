
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>

#include "gemmi/align.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/mmread.hpp"
#include "gemmi/model.hpp"

std::string extractFileName(const std::string& fullPath);
void convert(const std::string& inputFileName, std::ostream& output);

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> [<output_file>]" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];

    if (argc == 3) {
        std::ofstream outfile(argv[2], std::ios_base::app);
        if (!outfile.is_open()) {
            std::cerr << "Error: cannot open for output: " << argv[2] << std::endl;
            return 1;
        }
        convert(inputFileName, outfile);
    } else {
        convert(inputFileName, std::cout);
    }

    return 0;
}

std::string extractFileName(const std::string& fullPath)
{
    size_t lastSlash = fullPath.find_last_of("/\\");
    std::string filenameWithExt = fullPath.substr((lastSlash == std::string::npos) ? 0 : lastSlash + 1);
    size_t dotPos = filenameWithExt.find_last_of('.');
    return (dotPos == std::string::npos) ? filenameWithExt : filenameWithExt.substr(0, dotPos);
}

void convert(const std::string& inputFileName, std::ostream& output)
{
    std::unordered_set<char> AA_set { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y' };
    gemmi::Structure st = gemmi::read_structure_file(inputFileName);
    std::string filename = extractFileName(inputFileName);

    if (st.name == "XXXX" && filename.size() >= 4) {
        std::string id_name_file = filename.substr(0, 4);
        std::transform(id_name_file.begin(), id_name_file.end(), id_name_file.begin(), [](unsigned char c) { return std::toupper(c); });
        st.name = id_name_file;
    }

    gemmi::assign_label_seq_id(st, false);

    std::string method = "";
    if (!st.meta.experiments.empty()) {
        method = st.meta.experiments.back().method;
    }
    std::string resolution = std::to_string(st.resolution);

    if (!st.models.empty()) {
        gemmi::Model& mod = st.models.front();

        remove_alternative_conformations(mod);
        if (!mod.chains.empty()) {
            for (const gemmi::Chain& chain : mod.chains) {
                if (!chain.residues.size() || gemmi::check_polymer_type(chain.get_polymer()) != gemmi::PolymerType::PeptideL) {
                    continue;
                }
                std::string fasta;
                std::string out;
                bool firstPolymerResidue = true;
                gemmi::SeqId::OptionalNum oldLabel_seq;

                for (const gemmi::Residue& res : chain.residues) {
                    if (res.entity_type != gemmi::EntityType::Polymer)
                        continue;

                    if (firstPolymerResidue) {
                        oldLabel_seq = res.label_seq;
                        firstPolymerResidue = false;
                    }

                    int gap = 0;
                    gemmi::SeqId::OptionalNum labelSeqId = res.label_seq;
                    if (oldLabel_seq.has_value() && labelSeqId.has_value() && std::abs(*labelSeqId - *oldLabel_seq) > 1) {
                        gap = std::abs(*labelSeqId - *oldLabel_seq);
                    }

                    if (res.entity_type == gemmi::EntityType::Polymer) {
                        char c = std::toupper(gemmi::find_tabulated_residue(res.name).one_letter_code);
                        if (AA_set.find(c) == AA_set.end()) {
                            c = 'X';
                        }

                        if (gap) {
                            gemmi::Entity* ent = gemmi::find_entity(res.subchain, st.entities);
                            for (int currentId = 1; currentId < gap; ++currentId) {
                                if (ent) {
                                    out += std::tolower(gemmi::find_tabulated_residue(ent->full_sequence[*oldLabel_seq + currentId - 1]).one_letter_code);
                                } else {
                                    out += 'X';
                                }
                            }
                        }
                        out += c;
                        oldLabel_seq = labelSeqId;
                    }
                }
                int nb_X = std::count_if(out.begin(), out.end(), [](char ch) { return ch == 'X' || (ch >= 'a' && ch <= 'z'); });
                if (out.size() - nb_X > 5) {
                    fasta += '>';
                    fasta += st.name;
                    fasta += chain.name;
                    fasta += '\t';
                    fasta += method;
                    fasta += '\t';
                    fasta += resolution;
                    fasta += '\n';
                    fasta += out;

                    output << fasta << "\n";
                }
            }
        }
    }
}