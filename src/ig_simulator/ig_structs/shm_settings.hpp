#pragma once

#include "../include_me.hpp"

enum SHMAction { SubstitutionSHM, DeletionSHM, InsertionSHM };

class SHM {
    size_t position_;
    SHMAction action_;

    // if substitution
    char new_nucleotide_;

    // if deletion
    size_t num_deleted_;

    // if insertion
    string insertion_;

public:
    SHM(size_t position, SHMAction action) :
            position_(position),
            action_(action),
            new_nucleotide_('A'),
            num_deleted_(0),
            insertion_() { }

    bool IsSubstitution() const { return action_ == SubstitutionSHM; }

    bool IsDeletion() const { return action_ == DeletionSHM; }

    bool IsInsertion() const { return action_ == InsertionSHM; }

    char Substitution() const { return new_nucleotide_; }

    size_t NumDeleted() const { return num_deleted_; }

    string InsertedString() const { return insertion_; }
};

class SHMSettings {
    vector<SHM> mutations_;

public:
    void Add(SHM mutation) {
        mutations_.push_back(mutation);
    }

    size_t Size() { return mutations_.size(); }

    typedef vector<SHM>::iterator shm_iterator;

    shm_iterator begin() { return mutations_.begin(); }

    shm_iterator end() { return mutations_.end(); }
};