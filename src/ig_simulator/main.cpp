#include "gene_database.hpp"
#include "vdj_recombinator.hpp"

int main(int argc, char *argv[]) {
    // test input data
    string vgenes_filename = string(argv[1]);
    string dgenes_filename = string(argv[2]);
    string jgenes_filename = string(argv[3]);

    // HC DB
    HC_GenesDatabase hc_database;
    hc_database.AddGenesFromFile(variable_gene, vgenes_filename);
    hc_database.AddGenesFromFile(diversity_gene, dgenes_filename);
    hc_database.AddGenesFromFile(join_gene, jgenes_filename);
    hc_database.Print(cout);

    // IgSimulator pipeline
    HC_SimpleRecombinator hc_vdj_recombinator(hc_database, HC_SimpleRecombinationCreator(hc_database));
    for(size_t i = 0; i < 10; i++) {
        auto vdj = hc_vdj_recombinator.CreateRecombination();
        vdj.Print(cout);
        cout << "-----------------" << endl;
    }

    return 0;
}