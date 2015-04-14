#include "gene_database.hpp"
#include "repertoire.hpp"
#include "multiplicity_creator.hpp"
#include "config.hpp"
#include "base_repertoire_creation/vdj_recombinator.hpp"
#include "base_repertoire_creation/exonuclease_remover.hpp"
#include "base_repertoire_creation/p_nucleotides_creator.hpp"
#include "base_repertoire_creation/n_nucleotides_creator.hpp"
#include "base_repertoire_creation/cdr_labeler.hpp"
#include "mutated_repertoire_creation/shm_creator.hpp"

HC_GenesDatabase_Ptr CreateHCDatabase(HC_InputParams params) {
    HC_GenesDatabase_Ptr hc_database(new HC_GenesDatabase());
    hc_database->AddGenesFromFile(variable_gene, params.vgenes_fname);
    hc_database->AddGenesFromFile(diversity_gene, params.dgenes_fname);
    hc_database->AddGenesFromFile(join_gene, params.jgenes_fname);
    return hc_database;
}

HC_Repertoire_Ptr CreateBaseRepertoire(HC_InputParams params, HC_GenesDatabase_Ptr hc_database) {
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
    HC_Repertoire_Ptr base_repertoire(new HC_Repertoire());

    // base multiplicity creator
    double base_lambda = double(params.basic_repertoire_params.base_repertoire_size) /
            params.basic_repertoire_params.mutated_repertoire_size;
    HC_ExponentialMultiplicityCreator base_multiplicity_creator(base_lambda);

    // cdr labeling
    HC_CDRLabelingStrategy cdr_labeling_strategy;
    HC_CDRLabeler cdr_labeler(cdr_labeling_strategy);

    for(size_t i = 0; i < params.basic_repertoire_params.base_repertoire_size; i++) {
        auto vdj = hc_vdj_recombinator.CreateRecombination();
        vdj = remover.CreateRemovingSettings(vdj);
        vdj = p_creator.CreatePNucleotides(vdj);
        vdj = n_creator.CreateNNucleotides(vdj);
        cout << (*vdj);
        cout << "-----------------" << endl;

        HC_VariableRegionPtr ig_variable_region = HC_VariableRegionPtr(new HC_VariableRegion(vdj));
        ig_variable_region = cdr_labeler.LabelCDRs(ig_variable_region);

        size_t multiplicity = base_multiplicity_creator.AssignMultiplicity(ig_variable_region);
        base_repertoire->Add(HC_Cluster(ig_variable_region, multiplicity));
    }
    return base_repertoire;
}

HC_Repertoire_Ptr CreateMutatedRepertoire(HC_InputParams params, HC_Repertoire_Ptr base_repertoire) {
    HC_Repertoire_Ptr mutated_repertoire(new HC_Repertoire());

    // mutated multiplicity creator
    double mutated_lambda = double(base_repertoire->NumberAntibodies()) /
            double(params.basic_repertoire_params.final_repertoire_size);
    HC_ExponentialMultiplicityCreator mutated_multiplicity_creator(mutated_lambda);

    // shm creator
    HC_RgywWrcySHMStrategy shm_creation_strategy1(params.pattern_shm_params.min_number_pattern_shm,
                                                 params.pattern_shm_params.max_number_pattern_shm,
                                                 params.pattern_shm_params.substitution_propability);
    HC_CDRBasedRandomSHMStrategy shm_creation_strategy2(params.cdr_shm_params.min_number_mutations,
                                                        params.cdr_shm_params.max_number_mutations,
                                                        params.cdr_shm_params.mutation_in_fr_prop);
    HC_CompositeSHMCreationStrategy composite_shm_strategy(shm_creation_strategy1, shm_creation_strategy2);
    HC_SHMCreator shm_creator(composite_shm_strategy);

    for(auto it = base_repertoire->begin(); it != base_repertoire->end(); it++) {
        cout << "New base antibody, multiplicity: " << it->Multiplicity() << endl;
        for(size_t i = 0; i < it->Multiplicity(); i++) {
            auto variable_region_ptr = it->IgVariableRegion()->Clone();

            cout << "Old sequence: " << variable_region_ptr->Sequence() << endl;

            variable_region_ptr = shm_creator.CreateSHM(variable_region_ptr);
            size_t multiplicity = mutated_multiplicity_creator.AssignMultiplicity(variable_region_ptr);
            mutated_repertoire->Add(HC_Cluster(variable_region_ptr, multiplicity));

            cout << "New sequence: " << variable_region_ptr->Sequence() << endl;
            cout << "----" << endl;
        }
        cout << "-----------------------------" << endl;
        return mutated_repertoire;
    }
    return mutated_repertoire;
}

int CreateHCRepertoire(HC_InputParams params) {
    auto hc_database = CreateHCDatabase(params);
    HC_Repertoire_Ptr base_repertoire = CreateBaseRepertoire(params, hc_database);

    cout << "Base repertoire:" << endl;
    for(auto it = base_repertoire->begin(); it != base_repertoire->end(); it++) {
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;
        cout << it->IgVariableRegion()->GetCDRSettings();
    }

    HC_Repertoire_Ptr mutated_repertoire = CreateMutatedRepertoire(params, base_repertoire);

    cout << "Mutated repertoire:" << endl;
    for(auto it = mutated_repertoire->begin(); it != mutated_repertoire->end(); it++) {
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;
    }
    return 0;
}