#pragma once

#include "ig_structs/ig_structure_structs.hpp"

template<class GeneDB, class VDJ_Recombination_Ptr, class RecombinationCreator>
class VDJ_Recombinator {
    const GeneDB &gene_db_;
    RecombinationCreator recombination_creator_;

public:
    VDJ_Recombinator(const GeneDB &gene_db,
            RecombinationCreator recombination_creator) :
            gene_db_(gene_db),
            recombination_creator_(recombination_creator) { }

    VDJ_Recombination_Ptr CreateRecombination() {
        return recombination_creator_.CreateRecombination();
    }
};

// ----------------------------------------------------------------------------
//      Simple creator of VDJ recombination
//      - randomly chooses V, D, J genes from database
//      - probabilities from all V, D, and J are equal
// ----------------------------------------------------------------------------

// HC
class HC_SimpleRecombinationCreator {
    const HC_GenesDatabase &db_;

public:
    HC_SimpleRecombinationCreator(const HC_GenesDatabase &db) :
            db_(db) { }

    HC_VDJ_Recombination_Ptr CreateRecombination() {
        return HC_VDJ_Recombination_Ptr(
                new HC_VDJ_Recombination(db_,
                RandomIndex(db_.GenesNumber(variable_gene)),
                RandomIndex(db_.GenesNumber(diversity_gene)),
                RandomIndex(db_.GenesNumber(join_gene))));
    }
};

// LC
class LC_SimpleRecombinationCreator {
    const LC_GenesDatabase &db_;
    size_t seed_;

public:
    LC_SimpleRecombinationCreator(const LC_GenesDatabase &db) :
            db_(db) { }

    LC_VDJ_Recombination CreateRecombination() {
        return LC_VDJ_Recombination(db_,
                RandomIndex(db_.GenesNumber(variable_gene)),
                RandomIndex(db_.GenesNumber(join_gene)));
    }
};

// ----------------------------------------------------------------------------

typedef VDJ_Recombinator<HC_GenesDatabase, HC_VDJ_Recombination_Ptr, HC_SimpleRecombinationCreator>
        HC_SimpleRecombinator;
typedef VDJ_Recombinator<LC_GenesDatabase, LC_VDJ_Recombination_Ptr, LC_SimpleRecombinationCreator>
        LC_SimpleRecombinator;