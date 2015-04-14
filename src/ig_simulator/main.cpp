#include "config.hpp"
#include "hc_repertoire_launch.hpp"
#include "lc_repertoire_launch.hpp"

int main(int argc, char *argv[]) {
    HC_InputParams input_params;

    input_params.vgenes_fname =  string(argv[1]);
    input_params.dgenes_fname =  string(argv[2]);
    input_params.jgenes_fname =  string(argv[3]);

    input_params.basic_repertoire_params.base_repertoire_size = 10;
    input_params.basic_repertoire_params.mutated_repertoire_size = 100;
    input_params.basic_repertoire_params.final_repertoire_size = 750;

    input_params.pattern_shm_params = PatternSHMParams::CreateStandardParams();
    input_params.cdr_shm_params = CDR_SHMParams::CreateStandardParams();

    CreateHCRepertoire(input_params);
}