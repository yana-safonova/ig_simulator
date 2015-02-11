#pragma once

#include "gene_database.hpp"

class HC_VDJ_Recombination {
    const HC_GenesDatabase &db_;
    size_t vgene_index_;
    size_t dgene_index_;
    size_t jgene_index_;

public:
    HC_VDJ_Recombination(const HC_GenesDatabase &db,
        size_t vgene_index,
        size_t dgene_index,
        size_t jgene_index) :
            db_(db),
            vgene_index_(vgene_index),
            dgene_index_(dgene_index),
            jgene_index_(jgene_index) {
    }

    size_t VgeneIndex() { return vgene_index_; }

    size_t DgeneIndex() { return dgene_index_; }

    size_t JgeneIndex() { return jgene_index_; }

    void Print(ostream &out) {
        out << "Indices: " << vgene_index_ << " " << dgene_index_ << " " << jgene_index_ << endl;
        out << "Vgene. " << db_.GetByIndex(variable_gene, vgene_index_) << endl;
        out << "Dgene. " << db_.GetByIndex(diversity_gene, dgene_index_) << endl;
        out << "Jgene. " << db_.GetByIndex(join_gene, jgene_index_) << endl;
    }
};

// ----------------------------------------------------------------------------

class LC_VDJ_Recombination {
    const LC_GenesDatabase &db_;
    size_t vgene_index_;
    size_t jgene_index_;

public:
    LC_VDJ_Recombination(const LC_GenesDatabase &db,
            size_t vgene_index,
            size_t jgene_index) :
            db_(db),
            vgene_index_(vgene_index),
            jgene_index_(jgene_index) { }

    size_t VgeneIndex() { return vgene_index_; }

    size_t JgeneIndex() { return jgene_index_; }
};

// ----------------------------------------------------------------------------

template<class GeneDB, class VDJ_Recombination, class RecombinationCreator>
class VDJ_Recombinator {
    const GeneDB &gene_db_;
    RecombinationCreator recombination_creator_;

public:
    VDJ_Recombinator(const GeneDB &gene_db,
            RecombinationCreator recombination_creator) :
            gene_db_(gene_db),
            recombination_creator_(recombination_creator) { }

    VDJ_Recombination CreateRecombination() {
        return recombination_creator_.CreateRecombination();
    }
};

// ----------------------------------------------------------------------------
size_t RandomIndex(size_t to, size_t from = 0) {
    return rand() % to + from;
}

class HC_SimpleRecombinationCreator {
    const HC_GenesDatabase &db_;
    size_t seed_;

public:
    HC_SimpleRecombinationCreator(const HC_GenesDatabase &db, size_t seed = 1) :
            db_(db), seed_(seed) { }

    HC_VDJ_Recombination CreateRecombination() {
        size_t ind = RandomIndex(db_.GenesNumber(join_gene));
        return HC_VDJ_Recombination(db_,
                RandomIndex(db_.GenesNumber(variable_gene)),
                RandomIndex(db_.GenesNumber(diversity_gene)),
                RandomIndex(db_.GenesNumber(join_gene)));
    }
};

class LC_SimpleRecombinationCreator {
    const LC_GenesDatabase &db_;
    size_t seed_;

public:
    LC_SimpleRecombinationCreator(const LC_GenesDatabase &db, size_t seed = 1) :
            db_(db), seed_(seed) { }

    LC_VDJ_Recombination CreateRecombination() {
        LC_VDJ_Recombination(db_,
                RandomIndex(db_.GenesNumber(variable_gene)),
                RandomIndex(db_.GenesNumber(join_gene)));
    }
};

// ----------------------------------------------------------------------------

typedef VDJ_Recombinator<HC_GenesDatabase, HC_VDJ_Recombination, HC_SimpleRecombinationCreator> HC_SimpleRecombinator;
typedef VDJ_Recombinator<LC_GenesDatabase, LC_VDJ_Recombination, LC_SimpleRecombinationCreator> LC_SimpleRecombinator;