#pragma once

#include "ig_structure_structs.hpp"

template<class VDJ_Recombination_Ptr, class NInsertionStrategy, class NInsertionSettings>
class NNucleotidesCreator {
    NInsertionStrategy n_insertion_strategy_;

public:
    NNucleotidesCreator(NInsertionStrategy n_insertion_strategy) :
        n_insertion_strategy_(n_insertion_strategy) { }

    VDJ_Recombination_Ptr CreateNNucleotides(VDJ_Recombination_Ptr vdj_recombination) {
        NInsertionSettings settings = n_insertion_strategy_.CreateNInsertionSettings(vdj_recombination);
        vdj_recombination->AddNInsertionSettings(settings);
        return vdj_recombination;
    }
};

// ----------------------------------------------------------------------------

struct HC_SimpleNInsertionConstants {
    const static size_t vd_insertion_max_len = 10;
    const static size_t dj_insertion_max_len = 10;
};

class HC_SimpleNInsertionStrategy {
public:
    HC_NInsertionSettings CreateNInsertionSettings(HC_VDJ_Recombination_Ptr vdj_recombination) {
        size_t vd_insertion_len = RandomInt(HC_SimpleNInsertionConstants::vd_insertion_max_len);
        size_t dj_insertion_len = RandomInt(HC_SimpleNInsertionConstants::dj_insertion_max_len);
        return HC_NInsertionSettings(GetRandomSequence(vd_insertion_len), GetRandomSequence(dj_insertion_len));
    }
};

// create analogue for light chain

// ----------------------------------------------------------------------------

typedef NNucleotidesCreator<HC_VDJ_Recombination_Ptr, HC_SimpleNInsertionStrategy, HC_NInsertionSettings>
    HC_SimpleNNucleotidesCreator;