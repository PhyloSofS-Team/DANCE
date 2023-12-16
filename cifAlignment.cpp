// Standard libraries
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <valarray>
#include <filesystem>


// Definitions
#define GEMMI_WRITE_IMPLEMENTATION

// Gemmi libraries
#include "gemmi/fstream.hpp"
#include "gemmi/mmread.hpp"
#include "gemmi/model.hpp"
#include "gemmi/pirfasta.hpp"
#include "gemmi/polyheur.hpp"
#include "gemmi/align.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"


std::string read_pir_or_fasta(std::istream& is, std::string &id);
std::string AA_list = "ACDEFGHIKLMNPQRSTVWY";
// some functions and structure definitions/declarations

struct GapSection {
    int seqIndex;
    int start;    
    int end;     
}; 
std::vector<GapSection> removeIsolatedResidues(std::vector<std::string>& seqs, int continentSize, int isolationDistance);
gemmi::Atom reconstructO(const gemmi::Atom *, const  gemmi::Atom *, const  gemmi::Atom *);
void CheckSeqs(std::vector<std::string> names);
void getNextSeqPosition(int &seqPosition,  const std::string &seq,  const std::string &resname);

std::vector< std::vector<double> > alignAndComputeRMSD(
    std::vector<std::string>& seqs, 
    const std::vector<std::string>& names, 
    gemmi::Structure& tst,
    int len_seq, 
    bool use_weights, 
    bool Calpha
);
bool removeDuplicates(
    std::vector<std::string>& names,
    std::vector<std::string>& seqs,
    gemmi::Structure& tst,
    std::vector<std::vector<double>>& rmsd_mat,
    double similarity_threshold,
    bool output_removed,
    const std::string& refName,
    int len_seq
);

int getMaxScoreSequence(const std::vector<std::string>& names,
                        const std::vector<std::string>& seqs,
                        const std::string& consensus,
                        int len_seq,
                        int blosum62[22][22]);

typedef std::vector<std::vector<std::vector<double>>> Tensor3D; // for the raw coords;
typedef std::vector<std::vector<bool>> Tensor2D_bool; //for the mask;

extern int blosum62[22][22];
std::string blosumAA = "ARNDCQEGHILKMFPSTWYVX-";
int main(int argc, char** argv) {
    
    if (argc == 1) {
        std::cout << "No arguments provided. Use the -h flag for help." << std::endl;
        return 1; // Exit with error code
    }


    std::string inAln;
    std::string cifDir;
    float similarity_threshold = 0.1;
    int continentSize = 4;
    int isolationDistance = 15;
    int nb_ref = 1;
    bool Calpha = false,
        output_pdb = false,
        output_cif = false,
        output_aln = false,
        output_rmsd = false,
        output_raw_coords = false,
        output_removed = false,
        use_weights = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:d:cws:pfaroun:x:y:")) != -1) {
        switch (opt) {
            case 'h':
                std::cout << "Usage: " << argv[0] << " [options]\n";
                std::cout << "-i                Path to the aln file\n";
                std::cout << "-d                Path to the cif directory\n";
                std::cout << "-c                Centermass and alignment on Calpha only\n";
                std::cout << "-w                Enable weighted alignment\n";
                std::cout << "-s                Set RMSD similarity threshold for conformation removing (default: 0.1A)\n";
                std::cout << "-p                Enable output pdb file\n";
                std::cout << "-f                Enable output cif file\n";
                std::cout << "-a                Enable output aln file\n";
                std::cout << "-r                Enable output RMSD file option\n";
                std::cout << "-o                Enable output raw coords file option\n";
                std::cout << "-u                Enable output of removed sequences file information\n";
                std::cout << "-n                Set the number of references (default: 1)\n";
                std::cout << "-x                Set the continent size (strictly superior to)(default: 4)\n";
                std::cout << "-y                Set the isolation distance (superior or equal to)(default: 15)\n";
                exit(0);
                break;
            case 'i':
                if (optarg == nullptr) {
                    std::cerr << "Error: -i option requires an argument (path to the aln file).\n";
                    exit(1);
                }
                inAln = optarg;
                break;
            case 'd':
                if (optarg == nullptr) {
                    std::cerr << "Error: -d option requires an argument (path to the cif directory).\n";
                    exit(1);
                }
                cifDir = optarg;
                break;
            case 'c':
                Calpha = true;
                break;
            case 'w':
                use_weights = true;
                break;
            case 's':
                similarity_threshold = std::stof(optarg);
                break;
            case 'p':
                output_pdb = true;
                break;
            case 'f':
                output_cif = true;
                break;
            case 'a':
                output_aln = true;
                break;
            case 'r':
                output_rmsd = true;
                break;
            case 'o':
                output_raw_coords = true;
                break;
            case 'u':
                output_removed = true;
                break;
            case 'n':
                nb_ref = std::stoi(optarg); 
                break;
            case 'x': 
                continentSize = std::stoi(optarg);
                break;
            case 'y':  
                isolationDistance = std::stoi(optarg);
                break;
            default: 
                std::cerr << "Invalid option.\n";
                exit(1);
        }
        
    }

    std::string id;
    std::string seq;

    std::map <std::string, std::string> fasta_sequences;
    std::vector<std::string> names;
    std::vector<std::string> seqs;

    gemmi::Ifstream stream(inAln); // initial stream
    while(stream.ref()) {
        seq = read_pir_or_fasta(stream.ref(), id);
        fasta_sequences[id] = seq;
    }

    for (const auto& p : fasta_sequences ) {
        size_t pos = p.first.find('\t');  
        if (pos != std::string::npos) {
            names.push_back(p.first.substr(0, pos));
        } else {
            names.push_back(p.first);  
        }

        seqs.push_back(p.second);
    }      

    std::string refName=names[0];
    int len_seq = seq.size();
    
    std::cout<< "Phase 1 : collecting data"<<std::endl;
    // Phase 1 : collect backbone atoms, complete backbone-O and prune the MSA letters
    // first, check every chain for completeness and adjust the sequences accordingly
    gemmi::Structure tst;
    gemmi::Structure st;
    gemmi::Model &model = st.models[0];
    tst.models.reserve(seqs.size());
    //std::vector<std::string> methods;
    std::string cp_struct_file = "foo";
    for (size_t iSeq = 0; iSeq < names.size();) {
        std::string &seq = seqs[iSeq];
        const std::string &name = names[iSeq];

        gemmi::Chain tch("A");
        tch.residues.reserve(len_seq);

        std::string chain_name = name.substr(4);
        
        std::string name_transformed = name.substr(0,4);
        std::transform(name_transformed.begin(), name_transformed.end(), name_transformed.begin(),
            [](unsigned char c){ return std::tolower(c); });
        std::string struct_file = cifDir + name_transformed + "_final.cif";
        if (!(std::filesystem::exists(struct_file))){
            struct_file = cifDir+name.substr(0,4)+"_final.cif";
        }
        if (!(std::filesystem::exists(struct_file))){
            struct_file = cifDir+name_transformed + ".cif";
        }
        if (!(std::filesystem::exists(struct_file))){
            struct_file = cifDir+name.substr(0,4) + ".cif";
        }
        if (!(std::filesystem::exists(struct_file))){
            struct_file = cifDir+name.substr(0,4) + chain_name + ".cif";
        }
        if (!(std::filesystem::exists(struct_file))){
            std::string chain_name_transformed = chain_name;
            std::transform(chain_name_transformed.begin(), chain_name_transformed.end(), chain_name_transformed.begin(),
                [](unsigned char c){ return std::tolower(c); });
            struct_file = cifDir+name_transformed + chain_name_transformed + ".cif";
        }

        if (!(std::filesystem::exists(struct_file))){
            std::cout<<"No structure file found for "<<name<<std::endl;
            std::cout<<"Please check the input directory and the naming of the files"<<std::endl;
            exit(0);
        }


        

        if (struct_file != cp_struct_file){

            cp_struct_file= struct_file;
            std::cout<< "reading "<< struct_file<<std::endl;

            st = gemmi::read_structure_file(struct_file, gemmi::CoorFormat::Mmcif);
            gemmi::Model &model = st.models[0];
            remove_alternative_conformations(model);

        }

        gemmi::Model &model = st.models[0];
        gemmi::Chain* chain = model.find_chain(chain_name);
        std::cout<< "reading chain "<< chain_name<<std::endl;
        int seqPosition = -1;
        for(auto const &res : chain->residues) { // iterate over residues
            // this part checks the backbone and completes it if needed, the complete one is written to tch.residues
            if (res.entity_type == gemmi::EntityType::Polymer && AA_list.find(std::toupper(gemmi::find_tabulated_residue(res.name).one_letter_code)) != std::string::npos) {
                getNextSeqPosition(seqPosition, seq, res.name);

                gemmi::Residue tres;
                tres.atoms.reserve(8);
                tres.name = res.name;
                tres.entity_type = res.entity_type;
                tres.label_seq = res.label_seq;
                tres.seqid.num = seqPosition+1;

                const gemmi::Atom *atomN, *atomCa, *atomC, *atomO;

                atomN = res.get_n();
                atomCa = res.get_ca();
                atomC = res.get_c();
                if(atomN == NULL || atomCa == NULL || atomC == NULL) {

                    std::cout<< "backbone atoms missing"<<std::endl;
                    //remove current letter from the sequence in the MSA
                    seq[seqPosition] = '-';
                    continue; // backbone is incomplete
                }

                atomO = res.find_atom("O", '*', gemmi::El::O);
                if(atomO == NULL) { // no O in the backbone, try to reconstruct it...
                    const gemmi::Residue* nextRes = chain->next_residue(res);
                    if(nextRes==NULL) {

                        std::cout<< "next res missing"<<std::endl;

                        //remove current letter from the sequence in the MSA
                        seq[seqPosition] = '-';
                        continue; // no next residue
                    }

                    const gemmi::Atom *nextN = nextRes->get_n();
                    if(nextN==NULL) {

                        std::cout<< "next res N missing"<<std::endl;

                        //remove current letter from the sequence in the MSA
                        seq[seqPosition] = '-';
                        continue; // backbone is incomplete
                    }

                    gemmi::Atom atom = reconstructO(atomCa, atomC, nextN); // reconstruct O

                    std::cout<< "O reconstructed"<<std::endl;

                    tres.atoms.push_back(*atomN);
                    tres.atoms.push_back(*atomCa);
                    tres.atoms.push_back(*atomC);
                    tres.atoms.push_back(atom);
                } else {

                    tres.atoms.push_back(*atomN);
                    tres.atoms.push_back(*atomCa);
                    tres.atoms.push_back(*atomC);
                    tres.atoms.push_back(*atomO);
                }

                //backbone is complete!
                tch.residues.push_back(tres);
            } // end of the completion part, the MSA strings are adjusted accordingly!
        } // res iterator
        if (tch.get_polymer().length() >= 5.0 ) { // it's allright, save it in the multimodel tst
            gemmi::Model tmo(name);
            tmo.chains.push_back(tch);
            tst.models.push_back(tmo);
            //methods.push_back(st.meta.experiments.back().method);
            iSeq++;

        } else { //otherwise, remove the corresponding string from the MSA!

            std::cout<< "removing short sequence "<<names[iSeq]<<std::endl;
            names.erase(names.begin()+iSeq);
            seqs.erase(seqs.begin()+iSeq);

        }
    }

    

    // Phase 2 : Checking alignement:
    CheckSeqs(names); 
    std::cout<< "Phase 2 : Checking sequence alignement and chosing reference"<<std::endl;
    //Here we remove the isolated residues
    std::vector<GapSection> gapSections = removeIsolatedResidues(seqs, continentSize, isolationDistance);
    for (const auto& section : gapSections) {
        int seqIndex = section.seqIndex;
        int start = section.start;
        int end = section.end;

        gemmi::Model& model = tst.models[seqIndex]; 
        auto& residues = model.chains[0].residues; 

        auto it = residues.begin();
        while (it != residues.end()) {
            if (it->seqid.num.value >= start+1 && it->seqid.num.value < end+1) {
                std::cout << "Deleting residue at position " << it->seqid.num.value 
                        << " in sequence " << names[seqIndex] << std::endl;
                it = residues.erase(it);
            } else {
                ++it;
            }
        }
    }
    //Here we check if there are problematic columns (only X and -)
    for (int i=len_seq-1; i>=0; i--){
        bool nores = true;
        for (size_t j=0 ; j< names.size(); ++j) {
            if (seqs[j][i]!='X' && seqs[j][i] != '-'){
                nores=false;
                break;
            }
        }
        if (nores){
            std::cout<<"removing problematic column #"<<i<<std::endl;
            for(auto& r : seqs) r.erase(r.begin()+i);
            len_seq = len_seq-1;
            for (gemmi::Model &mod : tst.models){
                for(gemmi::Residue &res : mod.chains[0].residues){
                    if (res.seqid.num.value>i){res.seqid.num.value=res.seqid.num.value-1;}
                }
            }
        }
    }


std::string consensus = "";
for (int i = 0; i < len_seq; ++i) {
    std::map<char, int> freq_map;
    for (size_t seq_idx = 0; seq_idx < names.size(); ++seq_idx) {
        char current_aa = seqs[seq_idx][i];
        freq_map[current_aa]++;
    }

    char max_aa = seqs[0][i]; 
    int max_count = freq_map[max_aa];

    for (auto& entry : freq_map) {
        char aa = entry.first;
        int count = entry.second;

        if (count > max_count) {
            max_count = count;
            max_aa = aa;
        } 
        else if (count == max_count) {
            int aa_idx = blosumAA.find(aa);
            int max_aa_idx = blosumAA.find(max_aa);

            if (blosum62[aa_idx][aa_idx] < blosum62[max_aa_idx][max_aa_idx] && blosum62[aa_idx][aa_idx] > 0) {
                max_aa = aa;
            } 
            else if (blosum62[aa_idx][aa_idx] == blosum62[max_aa_idx][max_aa_idx]) {
                if (aa < max_aa) { //alphabetical order
                    max_aa = aa;
                }
            }
            else if (blosum62[max_aa_idx][max_aa_idx] == 0 && blosum62[aa_idx][aa_idx] > 0) {
                max_aa = aa;
            }
        }
        
    }

    consensus += max_aa;
}



    /* old method with blosum score
    std::vector<int> scores(names.size(), 0); // Create a score vector initialized with zeros.

    // For each sequence, calculate the BLOSUM62 substitution score with all other sequences.
    for (size_t seq_idx = 0; seq_idx < names.size(); ++seq_idx) {
        for (int i = 0; i < len_seq; ++i) { // Iterate over each amino acid in the sequence.
            char current_aa = seqs[seq_idx][i];
            
            int current_aa_idx = blosumAA.find(current_aa);
            
            for (size_t other_seq_idx = 0; other_seq_idx < names.size(); ++other_seq_idx) {
                if (seq_idx == other_seq_idx) continue; // Don't compare with itself.

                char other_aa = seqs[other_seq_idx][i];
                int other_aa_idx = blosumAA.find(other_aa);

                // Add the BLOSUM62 score to the current sequence score.
                scores[seq_idx] += blosum62[current_aa_idx][other_aa_idx];
            }
        }
    }

    // Find the sequence with the highest score.
    int max_score = scores[0];
    int maxarg = 0;
    for (size_t i = 1; i < names.size(); ++i) {
        if (scores[i] > max_score) {
            max_score = scores[i];
            maxarg = i;
        }
    }
*/
    /*//old method with score
    std::valarray<int> counterG(0,names.size());
   
    //find the reference
    for (int i=0 ; i<len_seq; ++i){ //iteration on AA
        std::valarray<int> counterv(0,names.size()); 
        int count = 0;
        for (size_t j=0 ; j< names.size(); ++j) { //iteration on column
            if (seqs[j][i]!='X' && seqs[j][i]!='-') {
                count = count+1;
                counterv[j]=1;
            }
        }
        if (count==1){ //penalty for the seq when residue is alone
            counterv=counterv*-names.size();
        }
        else{counterv=counterv*(count-1);}
        count=0;
        counterG =counterG + counterv;
    }

    int maxref = 0;
    int maxarg = 0;
    for (size_t i=0 ; i< names.size();i++){
        if (counterG[i] > maxref){
            maxref = counterG[i];
            maxarg = i;
        }
    }*/
    
    // let's have the referenceId at the first position in the vectors!


    // here copy 
    int maxarg = 0;
    decltype(alignAndComputeRMSD(seqs, names, tst, len_seq, use_weights, Calpha)) rmsd_mat;

    std::set<int> selectedRefs;  // To keep track of already selected references
    std::vector<int> referenceOrder; // To store the order of references
    // A mapping to keep track of the original positions
    std::vector<int> originalPositions(names.size());

    for(int ref = 0; ref < nb_ref; ref++) {
        
        // For the first reference, get the sequence with the maximum score and make it the reference
         if(ref == 0) {
            // Operations specific to the first iteration
            maxarg = getMaxScoreSequence(names, seqs, consensus, len_seq, blosum62);
            std::swap(names[0], names[maxarg]);
            std::swap(seqs[0], seqs[maxarg]);
            std::swap(tst.models[0], tst.models[maxarg]);
            std::cout << "Reference Id: " << names[0] << std::endl;

            // Print the consensus sequence
            std::cout << "Consensus Sequence: \n" << consensus << std::endl;
            std::cout << "Chosen Sequence (" << names[0] << "): \n" << seqs[0] << std::endl;
         }
        
                
        if(ref == 1) {
            
            originalPositions.clear();
            originalPositions.resize(names.size());
            for (int i = 0; i < names.size(); i++) {
                originalPositions[i] = i;
            }
            selectedRefs.insert(0);
            referenceOrder.push_back(0);
            
            std::cout << "Selecting " << nb_ref - 1 << " references..." << std::endl;

            for(int currentRef = 1; currentRef < nb_ref; currentRef++) {
                double highestValue = 0.0;  
                int maxarg = 0;
                for(int i = 1; i < names.size(); i++) {

                    if(selectedRefs.find(originalPositions[i]) != selectedRefs.end()) {
                        // Skip the sequence if it's already been selected as a reference
                        continue;
                    }

                    // Compute average RMSD to all prior references
                    double value = 0.0;
                    for(const int & refIdx : selectedRefs) {
 
                        value += rmsd_mat[refIdx][i];
                    }

                    value /= selectedRefs.size();

                    if(value > highestValue) {
                
                
                        highestValue = value;
                        maxarg = i;
                    }
                }
                selectedRefs.insert(originalPositions[maxarg]);
                referenceOrder.push_back(originalPositions[maxarg]);
                if (selectedRefs.size() == nb_ref) {
                    break;
                }
            }
        }
        if(ref > 0) {
            // For iterations after the first reference, move the previous reference back to its original position
            int prevRefOriginalPos = originalPositions[0];
            
            std::swap(names[0], names[prevRefOriginalPos]);       
            std::swap(seqs[0], seqs[prevRefOriginalPos]);      
            std::swap(tst.models[0], tst.models[prevRefOriginalPos]);    
            std::swap(originalPositions[0], originalPositions[prevRefOriginalPos]);

            // Now, swap the next reference to position 0
            maxarg = referenceOrder[ref];
            int nextRefOriginalPos = originalPositions[maxarg];
            std::swap(names[0], names[maxarg]);
            std::swap(seqs[0], seqs[maxarg]);    
            std::swap(tst.models[0], tst.models[maxarg]);
            std::swap(originalPositions[0], originalPositions[maxarg]);
            referenceOrder.push_back(nextRefOriginalPos);
            std::cout << "Reference Id for model n° " << ref << ": " << names[0] << std::endl;
    }


        // Phase 3 : Aligning and superposing
        std::cout<< "Phase 3 : Aligning and superposing"<<std::endl;
        
        if (!Calpha){
            gemmi::CenterOfMass centermass = gemmi::calculate_center_of_mass(tst.models[0].chains[0]);
            //computes the center of mass of the reference with real weights (but sometimes it is NaN for some reason, in this case we use the Calphas)
            if(centermass.get().x!=centermass.get().x || centermass.get().y!=centermass.get().y || centermass.get().z!=centermass.get().z){
                Calpha = true;
                std::cout<<"NaN Problem detected, Calpha set to true!"<<std::endl;
            }
            else{
                for (gemmi::Residue &res : tst.models[0].chains[0].residues){
                    for (gemmi::Atom &atom : res.atoms){
                        atom.pos = atom.pos - centermass.get();
                    }
                }
        }
        }

        
        if(Calpha){
            //recalculating the center of mass using the Calphas
            double sx=0.0,sy=0.0,sz=0.0;
            int cacount = 0;
            for(gemmi::Residue &res : tst.models[0].chains[0].residues){
                sx+=res.get_ca()->pos.x;
                sy+=res.get_ca()->pos.y;
                sz+=res.get_ca()->pos.z;
                cacount+=1;}
            sx=sx/cacount;
            sy=sy/cacount;
            sz=sz/cacount;
            gemmi::Position centermass = gemmi::Position {sx,sy,sz};
            for (gemmi::Residue &res : tst.models[0].chains[0].residues){
                    for (gemmi::Atom &atom : res.atoms){
                        atom.pos = atom.pos - centermass;
                    }
            }
        }
        
        std::cout<<"Centermass Computed"<<std::endl;
            
        rmsd_mat = alignAndComputeRMSD(seqs, names, tst, len_seq, use_weights, Calpha); //structural aln is done here

        
        if(ref==0){
            std::cout<<"Phase 4 : Checking for duplicates"<<std::endl;
            bool sequencesRemoved = removeDuplicates(names, seqs, tst, rmsd_mat, similarity_threshold, output_removed, refName, len_seq);
        
            int iterationCountRem = 1;

            while (sequencesRemoved && use_weights) {
                iterationCountRem++;
                std::cout<<"Realigning and removing duplicates, iteration "<<iterationCountRem<<std::endl;
                rmsd_mat = alignAndComputeRMSD(seqs, names, tst, len_seq, use_weights, Calpha); //structural aln is done here
                sequencesRemoved = removeDuplicates(names, seqs, tst, rmsd_mat, similarity_threshold, output_removed, refName, len_seq); //removes duplicates is done here
            }

            CheckSeqs(names); 
        }
       
        

        std::cout<< "Phase 5 : Profit!"<<std::endl;

        gemmi::remove_anisou(tst);

        if (output_pdb){
            gemmi::PdbWriteOptions opts;
            opts.ter_records = false;
            tst.raw_remarks.push_back("SEQLEN    " + std::to_string(len_seq));
            for (size_t iSeq=0 ; iSeq< names.size(); ++iSeq) {
                tst.raw_remarks.push_back(std::to_string(iSeq+1) +" : "+ names[iSeq]);
            }
            gemmi::Ofstream osPDB(refName+"_"+names[0]+"_mm.pdb", &std::cout);
            gemmi::write_pdb(tst, osPDB.ref(), opts);
            tst.raw_remarks.clear();

        }
        
        if (output_cif){
            std::ofstream os(refName+"_"+names[0]+".cif");
            gemmi::cif::write_cif_to_stream(os, gemmi::make_mmcif_document(tst));
        }
        
        if (output_rmsd){

            std::string rmsdFile = refName +"_"+names[0] + "_rmsd.txt";
            FILE *fptr = fopen(rmsdFile.c_str(),"w+");
            for (size_t iSeq=0 ; iSeq< names.size(); ++iSeq) {
                for (size_t jSeq=0 ; jSeq< names.size(); ++jSeq) {
                    fprintf(fptr,"%.3e ",rmsd_mat[iSeq][jSeq]);
                }
                fprintf(fptr,"\n");
            }
        }
        if (output_aln){

            std::string alnFile = refName +"_"+names[0] + "_aln.fa";
            FILE *alno = fopen(alnFile.c_str(),"w+");
            for (int i=0 ; i< names.size(); ++i) {
                fprintf(alno,">%s\n",names[i].c_str());
                fprintf(alno,"%s\n",seqs[i].c_str());

            }
        }

        if (output_raw_coords){
            if(!Calpha){
                Tensor3D at_rc(names.size(), std::vector<std::vector<double>>(len_seq*4, std::vector<double>(3)));
                Tensor2D_bool atom_mask(names.size(), std::vector<bool>(len_seq * 4, false));
                size_t model_idx = 0;  
                for (gemmi::Model &mod : tst.models) {
                    for (gemmi::Residue &res : mod.chains[0].residues) {
                        size_t atom_idx = 0;
                        for (gemmi::Atom &atom : res.atoms) {
                            at_rc[model_idx][(res.seqid.num.value-1)*4+atom_idx][0] = atom.pos.x;
                            at_rc[model_idx][(res.seqid.num.value-1)*4+atom_idx][1] = atom.pos.y;
                            at_rc[model_idx][(res.seqid.num.value-1)*4+atom_idx][2] = atom.pos.z;
                            atom_mask[model_idx][(res.seqid.num.value-1)*4+atom_idx] = true;
                            atom_idx++;
                        }
                    }
                model_idx++;  
                }
                void saveTensor(const Tensor3D &tensor, const std::string &filename);
                void saveTensorBool(const Tensor2D_bool &tensor, const std::string &filename);
                saveTensor(at_rc, refName+"_"+names[0]+"_raw_coords.bin");
                saveTensorBool(atom_mask, refName+"_"+names[0]+"_raw_coords_mask.bin");
            }

            else{
            
                Tensor3D at_rc(names.size(), std::vector<std::vector<double>>(len_seq, std::vector<double>(3)));
                Tensor2D_bool atom_mask(names.size(), std::vector<bool>(len_seq, false));
                size_t model_idx = 0;  
                for (gemmi::Model &mod : tst.models) {
                    for (gemmi::Residue &res : mod.chains[0].residues) {
                        const gemmi::Atom* ca_atom = res.get_ca();
                        if (ca_atom) { 
                            at_rc[model_idx][res.seqid.num.value-1][0] = ca_atom->pos.x;
                            at_rc[model_idx][res.seqid.num.value-1][1] = ca_atom->pos.y;
                            at_rc[model_idx][res.seqid.num.value-1][2] = ca_atom->pos.z;
                            atom_mask[model_idx][res.seqid.num.value-1] = true;
                        }
                    }
                    model_idx++;  
                }
                void saveTensor(const Tensor3D &tensor, const std::string &filename);
                void saveTensorBool(const Tensor2D_bool &tensor, const std::string &filename);
                saveTensor(at_rc, refName+"_"+names[0]+"_raw_coords_ca.bin");
                saveTensorBool(atom_mask, refName+"_"+names[0]+"_raw_coords_ca_mask.bin");
            }
        }
    }
    return 0;
}


///////////////////////////////////////////////////
// aux functions
///////////////////////////////////////////////////

std::vector<GapSection> removeIsolatedResidues(std::vector<std::string>& seqs, int continentSize, int isolationDistance) {
    std::vector<GapSection> modifiedSections; 
    for (size_t seqIndex = 0; seqIndex < seqs.size(); ++seqIndex) {
        std::string& seq = seqs[seqIndex];
        std::vector<std::pair<int, int>> sections;

        for (size_t i = 0; i < seq.size();) {
            if (seq[i] != '-' && seq[i] != 'X') {
                size_t start = i;
                while (i < seq.size() && seq[i] != '-' && seq[i] != 'X') {
                    i++;
                }
                sections.push_back({static_cast<int>(start), static_cast<int>(i)});
            } else {
                i++;
            }
        }

        for (const auto& section : sections) {
            int start = section.first;
            int end = section.second;
            bool isContinent = end - start > continentSize;
            bool isIsolated = true;

            if (!isContinent) {
                for (const auto& continent : sections) {
                    if (continent.second - continent.first > continentSize) {
                        int distanceStart = std::abs(start - continent.first);
                        int distanceEnd = std::abs(end - continent.second);
                        int distanceToContinentStart = std::abs(start - continent.second);
                        int distanceToContinentEnd = std::abs(end - continent.first);

                        if (distanceStart <= isolationDistance || distanceEnd <= isolationDistance ||
                            distanceToContinentStart <= isolationDistance || distanceToContinentEnd <= isolationDistance) {
                            isIsolated = false;
                            break;
                        }
                    }
                }

                if (isIsolated) {
                    std::fill(seq.begin() + start, seq.begin() + end, '-');
                    modifiedSections.push_back({static_cast<int>(seqIndex), start, end});
                }
            }
        }
    }

    return modifiedSections;
}

inline std::string read_pir_or_fasta(std::istream& is, std::string &id) {
  std::string sequence;
  std::string line;
  std::getline(is, line);
  if (line[0] != '>')
    std::cerr <<  "PIR/FASTA files start with '>'" << std::endl;
  const bool pir_format = (line.size() > 3 && line[3] == ';');
    id = line.substr(1);
    if(id[id.size()-1]== ';') {
        id = id.substr(0,id.size()-1);
    }
  bool ensure_pir_format = false;
  bool desc_line = pir_format;
  size_t discard_n = 0;
  while (is.peek() != '>' && std::getline(is, line)) {
    for (char c : line) {
      if (('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z') || c =='-') {
        if (islower(c) && (c!='-')){
                c='X';
            }
        sequence += c;
      } else if (!std::isspace(c) && c != '-') {
        if (desc_line) {
          ensure_pir_format = true;
          sequence.clear();
          break;
        }
        if (pir_format && c == '*')
          return sequence.substr(discard_n);
          //std::cerr <<  "unexpected character in sequence file: " << c << std::endl;
      }
    }
    if (desc_line)
      discard_n = sequence.size();
    desc_line = false;
  }
  if (ensure_pir_format)
      std::cerr <<  "PIR sequence must end with '*'" << std::endl;
  return sequence;
}

// reconstruct O based on coordinates of N, CA, C
gemmi::Atom  reconstructO(const gemmi::Atom *atomCa,const   gemmi::Atom *atomC, const  gemmi::Atom *atomN) {

    gemmi::Atom atomO = *atomC;
    atomO.element = gemmi::El::O;
    atomO.name = "O";
    atomO.serial = std::max(atomC->serial, atomCa->serial) + 1;

    // a typical C-O distance in 1.24 A
    double d = 1.24;
    gemmi::Vec3 CaC = (atomC->pos - atomCa->pos).normalized();
    gemmi::Vec3 NC = (atomN->pos - atomCa->pos).normalized();
    gemmi::Vec3 CO = (CaC + NC).normalized();
    atomO.pos = atomC->pos + gemmi::Position((d*CO));

    return atomO;
}



void getNextSeqPosition(int &seqPosition,  const std::string &seq,  const std::string &resname) {

    for(++seqPosition; seqPosition<seq.size() && (seq[seqPosition] == 'X' || seq[seqPosition] == '-'); ++seqPosition) {
        continue;
    }
    assert(seq[seqPosition] == std::toupper(gemmi::find_tabulated_residue(resname).one_letter_code));
}


void CheckSeqs(std::vector<std::string> names){
    if (names.size() <= 1){
        std::cout<<"Not enough conformations left, multimodel can't be built!"<< std::endl;
        exit(0);
    }
}
void saveTensor(const Tensor3D &tensor, const std::string &filename) {
    std::ofstream out(filename, std::ios::binary);
    if (out.is_open()) {
        size_t numModels = tensor.size();
        out.write(reinterpret_cast<const char*>(&numModels), sizeof(size_t));
        size_t numSeqs = (numModels > 0) ? tensor[0].size() : 0;
        out.write(reinterpret_cast<const char*>(&numSeqs), sizeof(size_t));
        size_t numCoords = (numSeqs > 0) ? tensor[0][0].size() : 0;
        out.write(reinterpret_cast<const char*>(&numCoords), sizeof(size_t));

        for (const auto &model : tensor) {
            for (const auto &seq : model) {
                out.write(reinterpret_cast<const char*>(seq.data()), seq.size() * sizeof(double));
            }
        }
    }
}

void saveTensorBool(const Tensor2D_bool &tensor, const std::string &filename) {
    std::ofstream out(filename, std::ios::binary);
    if (out.is_open()) {
        size_t numRows = tensor.size();
        out.write(reinterpret_cast<const char*>(&numRows), sizeof(size_t));
        size_t numCols = (numRows > 0) ? tensor[0].size() : 0;
        out.write(reinterpret_cast<const char*>(&numCols), sizeof(size_t));

        for (const auto &row : tensor) {
            for (bool value : row) {
                char byteValue = value ? 1 : 0; // Convert bool to char
                out.write(&byteValue, sizeof(char));
            }
        }
    }
}

std::vector< std::vector<double> > alignAndComputeRMSD(
    std::vector<std::string>& seqs, 
    const std::vector<std::string>& names, 
    gemmi::Structure& tst,
    int len_seq, 
    bool use_weights, 
    bool Calpha
){
    std::vector<double> coverage(len_seq, 0.0);

    if (use_weights) {
        std::cout<<"Computing Coverage"<<std::endl;
        

        for (size_t iSeq = 0; iSeq < names.size(); ++iSeq) {
            const std::string &seqi = seqs[iSeq];
            for (int k = 0; k < len_seq; k++) {
                if (seqi[k] != '-' && seqi[k] != 'X') {
                    coverage[k]++;
                }
            }
        }
        for (int k = 0; k < len_seq; k++) {
            coverage[k] = (coverage[k] / names.size());
        }
    }

    std::vector< std::vector<double> > rmsd_mat;
    rmsd_mat.resize(seqs.size(), std::vector<double>(seqs.size(), 0.0));

    for (size_t iSeq=0 ; iSeq< names.size(); ++iSeq) {
        //std::cout<<"Aligning " <<names[iSeq]<<std::endl;
        std::string &seqi = seqs[iSeq];

        for (size_t jSeq=iSeq ; jSeq< names.size(); ++jSeq) {
            std::string &seqj = seqs[jSeq];

            //collect alignment arrays
            std::vector<gemmi::Position> posi, posj;
            posi.reserve(4*len_seq);
            posj.reserve(4*len_seq);

            int indexi = -1;
            int indexj = -1;
            std::vector<double> subset_coverage;
            for(int k=0; k< len_seq; k++) {

                bool aligned = false;
                if(seqi[k] != '-' && seqi[k] != 'X') {
                    indexi ++;
                    aligned  = true;
                }

                if(seqj[k] != '-' && seqj[k] != 'X') {
                    indexj ++;
                } else {
                    aligned = false;
                }

                if(aligned) { // fill in the arrays

                    if (Calpha){
                        if (use_weights) {
                            subset_coverage.push_back(coverage[k]);
                        }
                        const gemmi::Atom* atom_i = tst.models[iSeq].chains[0].residues[indexi].get_ca();
                        posi.push_back(atom_i->pos);
                        const gemmi::Atom* atom_j = tst.models[jSeq].chains[0].residues[indexj].get_ca();
                        posj.push_back(atom_j->pos);
                    }
                    else{
                        if (use_weights) {
                            for (int atom_index = 0; atom_index < 4; ++atom_index) {
                                subset_coverage.push_back(coverage[k]);
                            }
                        }
                        for (const auto &a : tst.models[iSeq].chains[0].residues[indexi].atoms)  posi.push_back(a.pos);
                        for (const auto &a : tst.models[jSeq].chains[0].residues[indexj].atoms)  posj.push_back(a.pos);
                    }
                }
            }

            assert(posi.size() == posj.size());
            gemmi::SupResult result;
            if (use_weights){
                result = gemmi::superpose_positions(posi.data(), posj.data(), posi.size(), subset_coverage.data());
            }
            else{
                result = gemmi::superpose_positions(posi.data(), posj.data(), posi.size(), NULL);
            }

            if(iSeq == 0 && jSeq) { //superpose  jSeq on iSeq
                for (gemmi::Residue &res : tst.models[jSeq].chains[0].residues)
                  for (gemmi::Atom &atom : res.atoms)
                    atom.pos = gemmi::Position(result.transform.apply(atom.pos));
            }

            // and store RMSD
            rmsd_mat[iSeq][jSeq] = result.rmsd;
            rmsd_mat[jSeq][iSeq] = result.rmsd;

        }
    }
    return rmsd_mat;
}
bool removeDuplicates(
    std::vector<std::string>& names,
    std::vector<std::string>& seqs,
    gemmi::Structure& tst,
    std::vector<std::vector<double>>& rmsd_mat,
    double similarity_threshold,
    bool output_removed,
    const std::string& refName,
    int len_seq
){  
    static bool isFirstCall = true;
    std::set<int, std::greater<int>> pos_to_remove;

    FILE *remfile = nullptr;

    if (output_removed){
        std::string remfilen = refName+"_"+names[0] + "_duplicates.csv";
        if (!isFirstCall) {
            remfile = fopen(remfilen.c_str(),"a+");
        }
        else {
            remfile = fopen(remfilen.c_str(),"w+");
            fprintf(remfile,"removed,reference,rmsd\n");
        }
    }
    isFirstCall = false;
    for (size_t iSeq=0 ; iSeq< names.size(); ++iSeq){
        for (size_t jSeq=iSeq+1 ; jSeq< names.size(); ++jSeq) { 
            if(rmsd_mat[iSeq][jSeq]<similarity_threshold || (iSeq == 0 && rmsd_mat[iSeq][jSeq]!=rmsd_mat[iSeq][jSeq])){ //RMSD threshold is here ! We remove any seq that has nan rmsd with reference
                auto i_in_set = pos_to_remove.find(iSeq);
                auto j_in_set = pos_to_remove.find(jSeq);
                if (i_in_set !=pos_to_remove.end() || j_in_set != pos_to_remove.end()){
                    continue;
                    //if we already removed one or both member(s) of the selected pair we skip the iteration
                } 
                std::string &seqi = seqs[iSeq];
                std::string &seqj = seqs[jSeq];

                bool discard_i = true;
                bool discard_j = true; //standard state : both seqs can be discarded
                if (!(iSeq == 0 && rmsd_mat[iSeq][jSeq]!=rmsd_mat[iSeq][jSeq])){ //we go straight to removing j if it has nan with ref
                    for (int iAA=0 ; iAA<len_seq; ++iAA){
                        if ( (seqi[iAA]==seqj[iAA]) || ((seqi[iAA] == '-' || seqi[iAA] == 'X') && (seqj[iAA] == '-' || seqj[iAA] == 'X'))) {
                            //we consider the '-' and 'X' identical here
                            continue;
                        }
                        else if ((seqi[iAA] != '-' && seqi[iAA] != 'X') && (seqj[iAA] != '-' && seqj[iAA] != 'X')){ 
                            //AA are differents but not gap or X -> there is sequence variability so we keep them both
    
                            discard_i = false;
                            discard_j = false;
                            break;
                        }
                        else if ((seqi[iAA] != '-' && seqi[iAA] != 'X') && (seqj[iAA] == '-' || seqj[iAA] == 'X')){
                            //i has AA where j has not, we cannot discard it
                            discard_i = false;
                            continue;
                        }
                        else if ((seqj[iAA] != '-' && seqj[iAA] != 'X') && (seqi[iAA] == '-' || seqi[iAA] == 'X')){
                            //j has AA where i has not, we cannot discard it
                            discard_j = false;
                            continue;
                        }
                        else{
                            std::cout<<"Should not happen!!!"<<std::endl;
                            std::cout<<seqi[iAA]<<seqj[iAA]<<std::endl;
                        } 
    
                        
                    }
                }
                if (discard_j){
                    //we discard the last in alphabetical order if it can be
                    pos_to_remove.insert(jSeq);
                    std::cout<<"conformation n° "<<jSeq<<" ("<<names[jSeq]<<") is identical to - or included in n° "<<iSeq<<" ("<<names[iSeq]<<")"<<std::endl;
                    if (remfile){
                        fprintf(remfile,"%s,%s,%.3e\n",names[jSeq].c_str(),names[iSeq].c_str(),rmsd_mat[iSeq][jSeq]);
                    }
                }
                else if (discard_i && iSeq!=0){
                    pos_to_remove.insert(iSeq);
                    std::cout<<"conformation n° "<<iSeq<<" ("<<names[iSeq]<<") is identical to - or included in n° "<<jSeq<<" ("<<names[jSeq]<<")"<<std::endl;
                    if (remfile){
                        fprintf(remfile,"%s,%s,%.3e\n",names[iSeq].c_str(),names[jSeq].c_str(),rmsd_mat[iSeq][jSeq]);
                    }
                }
            }

        }
    }

    if (remfile) {
        fclose(remfile);
    }

    for (int i : pos_to_remove){
        std::cout<<"removing n° "<<i<<" ("<<names[i]<<")"<<std::endl;
        tst.remove_model(names[i]);
        names.erase(names.begin()+i);
        seqs.erase(seqs.begin()+i);
        //methods.erase(methods.begin()+i);
        for(auto& r : rmsd_mat) r.erase(r.begin()+i);
        rmsd_mat.erase(rmsd_mat.begin()+i);
    }
    return !pos_to_remove.empty();
}

int getMaxScoreSequence(const std::vector<std::string>& names,
                        const std::vector<std::string>& seqs,
                        const std::string& consensus,
                        int len_seq,
                        int blosum62[22][22]) {
    int max_score = INT_MIN;
    int maxarg = -1;

    for (size_t seq_idx = 0; seq_idx < names.size(); ++seq_idx) {
        int current_score = 0;
        for (int i = 0; i < len_seq; ++i) {
            char consensus_aa = consensus[i];
            char current_aa = seqs[seq_idx][i];

            int consensus_aa_idx = blosumAA.find(consensus_aa); // Ensure blosumAA is available within this scope
            int current_aa_idx = blosumAA.find(current_aa);

            current_score += blosum62[consensus_aa_idx][current_aa_idx];
        }

        if (current_score > max_score) {
            max_score = current_score;
            maxarg = seq_idx;
        }
    }
    return maxarg;
}
int blosum62[22][22] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   X   -
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5, -5}, //A
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5, -5}, //R
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5, -5}, //N
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5, -5}, //D
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5, -5}, //C
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5, -5}, //Q
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5, -5}, //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5, -5}, //G
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5, -5}, //H
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5, -5}, //I
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5, -5}, //L
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5, -5}, //K
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5, -5}, //M
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5, -5}, //F
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5, -5}, //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5, -5}, //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5, -5}, //T
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5, -5}, //W
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5, -5}, //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5, -5}, //V
{-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  0,  0}, //X
{-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  0,  0}  //-
};

