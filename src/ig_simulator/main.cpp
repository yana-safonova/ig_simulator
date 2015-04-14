#include "config.hpp"
#include "hc_repertoire_launch.hpp"
#include "lc_repertoire_launch.hpp"


/*
 * ./ig_simulator HC base_rep_size mutated_rep_size final_rep_size Vgene.fa Dgene.fa Jgene.fa
 * ./ig_simulator LC base_rep_size mutated_rep_size final_rep_size Vgene.fa Jgene.fa
 */

void HCUsage() {
    cout << "Usage for simulation heavy chain repertoire:" << endl;
    cout << "./ig_simulator HC base_rep_size mutated_rep_size final_rep_size Vgene.fa Dgene.fa Jgene.fa" << endl;
}

void LCUsage() {
    cout << "Usage for simulation light chain repertoire:" << endl;
    cout << "./ig_simulator LC base_rep_size mutated_rep_size final_rep_size Vgene.fa Jgene.fa" << endl;
}

void Usage() {
    HCUsage();
    LCUsage();
}

enum ChainType {Unknown_chain, Heavy_chain, Light_chain};

ChainType ParseChainType(string arg) {
    if(arg == "HC")
        return Heavy_chain;
    else if(arg == "LC")
        return Light_chain;
    else
        return Unknown_chain;
}

HC_InputParams ParseHCInputParams(int argc, char* argv[]) {
    HC_InputParams input_params;
    input_params.basic_repertoire_params.base_repertoire_size = StringToType<int>(string(argv[2]));
    input_params.basic_repertoire_params.mutated_repertoire_size = StringToType<int>(string(argv[3]));
    input_params.basic_repertoire_params.final_repertoire_size = StringToType<int>(string(argv[4]));
    input_params.vgenes_fname =  string(argv[5]);
    input_params.dgenes_fname =  string(argv[6]);
    input_params.jgenes_fname =  string(argv[7]);
    input_params.pattern_shm_params = PatternSHMParams::CreateStandardParams();
    input_params.cdr_shm_params = CDR_SHMParams::CreateStandardParams();
    return input_params;
}

LC_InputParams ParseLCInputParams(int argc, char* argv[]) {
    LC_InputParams input_params;
    input_params.basic_repertoire_params.base_repertoire_size = StringToType<int>(string(argv[2]));
    input_params.basic_repertoire_params.mutated_repertoire_size = StringToType<int>(string(argv[3]));
    input_params.basic_repertoire_params.final_repertoire_size = StringToType<int>(string(argv[4]));
    input_params.vgenes_fname =  string(argv[5]);
    input_params.jgenes_fname =  string(argv[6]);
    input_params.pattern_shm_params = PatternSHMParams::CreateStandardParams();
    input_params.cdr_shm_params = CDR_SHMParams::CreateStandardParams();
    return input_params;
}

int main(int argc, char *argv[]) {
    if(argc == 1) {
        cout << "ERROR: Invalid number of input parameters" << endl;
        Usage();
        return 1;
    }

    string chain_type_arg = string(argv[1]);
    ChainType chain_type = ParseChainType(chain_type_arg);
    if(chain_type == Unknown_chain) {
        cout << "ERROR: Invalid chain type: " << chain_type_arg << endl;
        Usage();
        return 1;
    }

    if(chain_type == Heavy_chain) {
        if(argc != 8) {
            cout << "ERROR: Invalid number of input parameters" << endl;
            HCUsage();
            return 1;
        }
        HC_InputParams input_params = ParseHCInputParams(argc, argv);
        CreateHCRepertoire(input_params);
    }
    else {
        if(argc != 7) {
            cout << "ERROR: Invalid number of input parameters" << endl;
            LCUsage();
            return 1;
        }
        LC_InputParams input_params = ParseLCInputParams(argc, argv);
        CreateLCRepertoire(input_params);
    }
}