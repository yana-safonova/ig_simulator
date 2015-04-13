#include "hc_repertoire_launch.hpp"

int main(int argc, char *argv[]) {
    HC_InputParams input_params;
    input_params.vgenes_fname =  string(argv[1]);
    input_params.dgenes_fname =  string(argv[2]);
    input_params.jgenes_fname =  string(argv[3]);

    input_params.base_repertoire_size = 10;
    input_params.mutated_repertoire_size = 100;
    input_params.final_repertoire_size = 750;

    CreateHCRepertoire(input_params);
}