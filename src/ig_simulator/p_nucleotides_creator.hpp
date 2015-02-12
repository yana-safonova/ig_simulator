#pragma once

#include "ig_structure_structs.hpp"

template<class VDJ_Recombination_Ptr, class PInsertionStrategy, class PInsertionSettings>
class PNucleotidesCreator {
    PInsertionStrategy p_insertion_strategy_;

public:
    PNucleotidesCreator(PInsertionStrategy p_insertion_strategy) :
            p_insertion_strategy_(p_insertion_strategy) { }

    VDJ_Recombination_Ptr CreatePNucleotides(VDJ_Recombination_Ptr vdj_recombination) {
        PInsertionSettings settings = p_insertion_strategy_.CreatePInsertionSettings(vdj_recombination);
        vdj_recombination->AddPInsertionSettings(settings);
        return vdj_recombination;
    }
};

// ----------------------------------------------------------------------------

struct HC_SimplePInsertionConstants {
    const static size_t v_end_max = 3;
    const static size_t d_start_max = 3;
    const static size_t d_end_max = 3;
    const static size_t j_start_max = 3;
};

class HC_SimplePInsertionStrategy {
public:
    HC_PInsertionSettings CreatePInsertionSettings(HC_VDJ_Recombination_Ptr recombination) {
        return HC_PInsertionSettings(
                RandomInt(HC_SimplePInsertionConstants::v_end_max),
                RandomInt(HC_SimplePInsertionConstants::d_start_max),
                RandomInt(HC_SimplePInsertionConstants::d_end_max),
                RandomInt(HC_SimplePInsertionConstants::j_start_max)
        );
    }
};

// add analogue for light chain

// ----------------------------------------------------------------------------

typedef PNucleotidesCreator<HC_VDJ_Recombination_Ptr, HC_SimplePInsertionStrategy, HC_PInsertionSettings>
    HC_SimplePNucleotidesCreator;