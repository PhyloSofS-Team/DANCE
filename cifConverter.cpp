
#include "libs/cParameters.hpp"
#include "libs/cTimer.hpp"
#include "libs/cInfo.hpp"
#include <string>
#include <fstream>

#include  <sstream>
#include <stdlib.h>
#include <iostream>

#include "gemmi/mmcif.hpp"
#include "gemmi/mmread.hpp"
#include "gemmi/model.hpp"
#include "gemmi/align.hpp"

void convert(cParameters &params, cInfo &info);
void writeFasta(cParameters &params, cInfo &info, const std::string &str);
void addFasta(cParameters &params, cInfo &info, const std::string &str);

/*********************************************************/
/*******here goes function and struct declarations********/
/*********************************************************/

void printHeader(int headerWidth = 75) {
	
	std::string str("\n cifConverter : \n\
e-mail: xxx"
					);

	std::istringstream iss(str);
	
	printf("*%.*s*\n",headerWidth,"******************************************************************************************************************************************************************");
	std::string line;
	while (std::getline(iss, line)) {
		int i1 = (headerWidth - line.size())/2;
		int i2=i1;
		if (i2*2+line.size() != headerWidth) i2++;
		printf("*%.*s",i1,"------------------------------------------------------------------------------------");
		printf("%s",line.c_str());
		printf("%.*s*\n",i2,"------------------------------------------------------------------------------------");
	}
	printf("*%.*s*\n",headerWidth,"------------------------------------------------------------------------------------------------------------------------------------------------------------------");
	printf("*%.*s*\n",headerWidth,"******************************************************************************************************************************************************************");
	
}

int main(int argc, char** argv) {

	// timer
	cSequentialTimer timer;


	// some output
	int len = 70;
	printHeader(len);

	cInfo       info(len);

	cParameters	params(&info);


	// parse the arguments
	params.parseInputArguments(argc, argv);
	info.saveTimeInterval(timer.getTimeInterval());

    convert(params, info);
    info.saveTimeInterval(timer.getTimeInterval());

    info.setTotalTime(timer.getTime());
    info.printTimings();

	return 0;
	
}

std::string extractFileName(const std::string& fullPath) {
    // Extract filename with extension from path
    size_t lastSlash = fullPath.find_last_of("/\\");
    std::string filenameWithExt = fullPath.substr((lastSlash == std::string::npos) ? 0 : lastSlash + 1);
    
    // Remove extension from filename
    size_t dotPos = filenameWithExt.find_last_of('.');
    return (dotPos == std::string::npos) ? filenameWithExt : filenameWithExt.substr(0, dotPos);
}

void convert(cParameters &params, cInfo &info) {

    std::ofstream outfile(params.outputFileName, std::ios_base::app);
    if (!outfile.is_open()) {
        info.printInfo("Error: cannot open for output", params.outputFileName);
        return; 
    }
    info.printInfo("Converting the files");
    info.printInfo("Input name", params.inputFileName);
    info.printInfo("Output name", params.outputFileName);
    std::string AA_list = "ACDEFGHIKLMNPQRSTVWXY";
    gemmi::Structure st = gemmi::read_structure_file(params.inputFileName);
    std::string filename = extractFileName(params.inputFileName);

    if (st.name == "XXXX" && filename.size() >= 4) {
        std::string id_name_file = filename.substr(0,4);
        std::transform(id_name_file.begin(), id_name_file.end(), id_name_file.begin(), [](unsigned char c){ return std::toupper(c); });
        st.name = id_name_file;
    }    
    
//    gemmi::setup_entities(st);

    // Uses sequence alignment (model to SEQRES) to assign label_seq.
    // force: assign label_seq even if full sequence is not known (assumes no gaps)
    gemmi::assign_label_seq_id(st, false);

    std::string method = "";
    if(st.meta.experiments.size())
        method = st.meta.experiments.back().method; // may not exist for vanilla PDB!
    std::string resolution = std::to_string(st.resolution); // may be 0 for vanilla PDB!

    // FIXME: only the first model is taken for multi-model PDBs (NMR, for example)
    if (st.models.size() > 0) {
        gemmi::Model &mod = st.models[0];

        // FIXME: we only treat the first conformation
        remove_alternative_conformations(mod);
        if (mod.chains.size()) {

            for (const gemmi::Chain &chain: mod.chains) {

                std::string fasta; // Fasta Header  ">SeqID
                std::string out; // output
                gemmi::SeqId firstResId;
                gemmi::SeqId::OptionalNum firstLabelResId;

                gemmi::SeqId oldResId;
                gemmi::SeqId::OptionalNum oldLabel_seq;

                if (chain.residues.size() && gemmi::check_polymer_type(chain.get_polymer()) == gemmi::PolymerType::PeptideL) {



                    bool firstPolymerResidue = true;
                    for (const gemmi::Residue &res: chain.residues) {

                        if (res.entity_type != gemmi::EntityType::Polymer) continue;

                        if (firstPolymerResidue) {
                            firstResId = res.seqid; // this still can be water
                            firstLabelResId = res.label_seq; // mmCIF _atom_site.label_seq_id

                            // the actual sequence
                            oldResId = firstResId;
                            oldLabel_seq = firstLabelResId;

                            firstPolymerResidue = false;
                        }

                        // res.label_seq is the result of assign_label_seq_id() OR mmCIF _atom_site.label_seq_id
                        int gap = 0;

                        gemmi::SeqId resId = res.seqid;
                        gemmi::SeqId::OptionalNum labelSeqId = res.label_seq;

                        if (std::abs(*labelSeqId - *oldLabel_seq) > 1) {
                            gap = std::abs(*labelSeqId - *oldLabel_seq);
                        }

                        if (res.entity_type == gemmi::EntityType::Polymer) {
                            oldResId = resId;

                            char c = std::toupper(gemmi::find_tabulated_residue(res.name).one_letter_code);

                            if (AA_list.find(c)==std::string::npos){
                                c = 'X';
                            }

                            if (gap) { // if the sequence is known, we put the sequence, otherwise just 'X'

                                // the case of PDB or PDB-Redo
                                gemmi::Entity *ent = gemmi::find_entity(res.subchain, st.entities); // entity of the residue subchain

                                for (int currentId = 1; currentId < gap; ++currentId) {
                                    if(ent) {
                                        out += std::tolower(gemmi::find_tabulated_residue(ent->full_sequence[oldLabel_seq.value+currentId-1]).one_letter_code);
                                    } else
                                        out += 'X';
                                }
                            }
                            out += c;
                            
                            oldLabel_seq = labelSeqId;

                        }

                    }

                    // FIXME :the condition has to be changed for small letters
                    int nb_X = 0;
                    for (char &ch:out)
                        if(ch=='X' || (ch>='a' && ch <= 'z')) {nb_X+=1;}

                    if(out.size()-nb_X>5){
                        fasta += '>';
                        fasta += st.name;
                        fasta += chain.name;
                        fasta += '\t';
                        fasta += method;
                        fasta += '\t';
                        fasta += resolution;
                        fasta += '\n';
                        fasta += out;
        //                writeFasta(params, info, fasta);
                        outfile << fasta << "\n";
                    }
                    else{info.printInfo("INFO : sequence removed - less than 5 AA in ", params.inputFileName);}
                }

                else {
                        info.printInfo("ERROR : zero-length chain or non peptide type in ", params.inputFileName);
                    }}
            }else {
                info.printInfo("ERROR : No chains in a model in ", params.inputFileName);
                return;
            }

    }else {
        info.printInfo("ERROR : No models in ", params.inputFileName);
        return;
    }
    outfile.close();
}

void writeFasta(cParameters &params, cInfo &info, const std::string &str) {

    FILE * pFile;
    pFile = fopen (params.outputFileName.c_str(),"w");

    if (pFile==NULL) {
        info.printInfo("Error : cannot open for output", params.outputFileName);
        return;
    }

    info.printInfo("Output filename", params.outputFileName);

    fprintf(pFile,"%s\n",str.c_str());

    fclose (pFile);
}
