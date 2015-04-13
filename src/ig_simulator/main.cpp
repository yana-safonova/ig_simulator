#include "gene_database.hpp"
#include "vdj_recombinator.hpp"
#include "exonuclease_remover.hpp"
#include "p_nucleotides_creator.hpp"
#include "n_nucleotides_creator.hpp"
#include "repertoire.hpp"
#include "multiplicity_creator.hpp"
#include "cdr_labeler.hpp"
//#include "shm_creator.hpp"

int main(int argc, char *argv[]) {
    // test input data
    string vgenes_filename = string(argv[1]);
    string dgenes_filename = string(argv[2]);
    string jgenes_filename = string(argv[3]);

    size_t base_repertoire_size = 10;
    size_t mutated_repertoire_size = 100;
    size_t final_repertoire_size = 750;

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

    // cdr labeling
    HC_CDRLabelingStrategy cdr_labeling_strategy;
    HC_CDRLabeler cdr_labeler(cdr_labeling_strategy);

    for(size_t i = 0; i < base_repertoire_size; i++) {
        auto vdj = hc_vdj_recombinator.CreateRecombination();
        vdj = remover.CreateRemovingSettings(vdj);
        vdj = p_creator.CreatePNucleotides(vdj);
        vdj = n_creator.CreateNNucleotides(vdj);
        cout << (*vdj);
        cout << "-----------------" << endl;

        HC_VariableRegionPtr ig_variable_region = HC_VariableRegionPtr(new HC_VariableRegion(vdj));

        ig_variable_region = cdr_labeler.LabelCDRs(ig_variable_region);

        size_t multiplicity = base_multiplicity_creator.AssignMultiplicity(ig_variable_region);
        base_repertoire.Add(HC_Cluster(ig_variable_region, multiplicity));
    }

    mutated_repertoire_size = 0;
    for(auto it = base_repertoire.begin(); it != base_repertoire.end(); it++)
        mutated_repertoire_size += it->Multiplicity();

    cout << "Base repertoire:" << endl;
    for(auto it = base_repertoire.begin(); it != base_repertoire.end(); it++) {
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;
        cout << it->IgVariableRegion()->GetCDRSettings();
    }

    // mutated repertoire
    HC_Repertoire mutated_repertoire;

    // mutated multiplicity creator
    double mutated_lambda = double(mutated_repertoire_size) / double(final_repertoire_size);
    HC_ExponentialMultiplicityCreator mutated_multiplicity_creator(mutated_lambda);

    //size_t min_number_mutations = 2;
    //size_t max_number_mutations = 10;
    //HC_RgywWrcySHMStrategy shm_creation_strategy(min_number_mutations, max_number_mutations);
    //HC_SHMCreator shm_creator(shm_creation_strategy);

    for(auto it = base_repertoire.begin(); it != base_repertoire.end(); it++) {
        cout << "New base antibody, multiplicity: " << it->Multiplicity() << endl;
        for(size_t i = 0; i < it->Multiplicity(); i++) {
            auto variable_region_ptr = it->IgVariableRegion()->Clone();

            //variable_region_ptr = shm_creator.CreateSHM(variable_region_ptr);
            size_t multiplicity = mutated_multiplicity_creator.AssignMultiplicity(variable_region_ptr);
            mutated_repertoire.Add(HC_Cluster(variable_region_ptr, multiplicity));
            cout << "----" << endl;
        }
        cout << "-----------------------------" << endl;
        return 0;
    }


    return 0;
}