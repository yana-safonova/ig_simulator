#include "gene_database.hpp"
#include "vdj_recombinator.hpp"
#include "exonuclease_remover.hpp"
#include "p_nucleotides_creator.hpp"
#include "n_nucleotides_creator.hpp"
#include "repertoire.hpp"
#include "multiplicity_creator.hpp"

int main(int argc, char *argv[]) {
    // test input data
    string vgenes_filename = string(argv[1]);
    string dgenes_filename = string(argv[2]);
    string jgenes_filename = string(argv[3]);

    size_t base_repertoire_size = 10;
    size_t mutated_repertoire_size = 100;

    // HC DB
    HC_GenesDatabase hc_database;
    hc_database.AddGenesFromFile(variable_gene, vgenes_filename);
    hc_database.AddGenesFromFile(diversity_gene, dgenes_filename);
    hc_database.AddGenesFromFile(join_gene, jgenes_filename);
    //hc_database.Print(cout);

    // IgSimulator pipeline
    HC_SimpleRecombinator hc_vdj_recombinator(hc_database, HC_SimpleRecombinationCreator(hc_database));

    // endonuclease remover
    HC_SimpleRemovingStrategy removing_strategy;
    HC_SimpleExonucleaseRemover remover(removing_strategy);

    // p nucleotides insertion
    HC_SimplePInsertionStrategy p_nucls_strategy;
    HC_SimplePNucleotidesCreator p_creator(p_nucls_strategy);

    // n nucleotides insertion
    HC_SimpleNInsertionStrategy n_nucls_strategy;
    HC_SimpleNNucleotidesCreator n_creator(n_nucls_strategy);

    // repertoire
    HC_Repertoire base_repertoire;

    // base multiplicity creator
    double base_lambda = double(base_repertoire_size) / mutated_repertoire_size;
    HC_ExponentialMultiplicityCreator base_multiplicity_creator(base_lambda);

    for(size_t i = 0; i < base_repertoire_size; i++) {
        auto vdj = hc_vdj_recombinator.CreateRecombination();
        vdj = remover.CreateRemovingSettings(vdj);
        vdj = p_creator.CreatePNucleotides(vdj);
        vdj = n_creator.CreateNNucleotides(vdj);
        cout << (*vdj);
        cout << "-----------------" << endl;

        HC_VariableRegionPtr ig_variable_region = HC_VariableRegionPtr(new HC_VariableRegion(vdj));
        size_t multiplicity = base_multiplicity_creator.AssignMultiplicity(ig_variable_region);
        base_repertoire.Add(HC_Cluster(ig_variable_region, multiplicity));
    }

    for(auto it = base_repertoire.begin(); it != base_repertoire.end(); it++)
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;

    return 0;
}