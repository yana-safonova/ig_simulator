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
            action_(action) { }

    size_t Position() { return position_; }

    bool IsSubstitution() const { return action_ == SubstitutionSHM; }

    bool IsDeletion() const { return action_ == DeletionSHM; }

    bool IsInsertion() const { return action_ == InsertionSHM; }

    // ----------------------

    void SetSubstitution(char new_nucleotide) {
        assert(action_ == SubstitutionSHM);
        new_nucleotide_ = new_nucleotide;
    }

    char GetSubstitution() const {
        assert(action_ == SubstitutionSHM);
        return new_nucleotide_;
    }

    // ----------------------

    void SetNumDeleted(size_t num_deleted) {
        assert(action_ == DeletionSHM);
        num_deleted_ = num_deleted;
    }

    size_t GetNumDeleted() const {
        assert(action_ == DeletionSHM);
        return num_deleted_;
    }

    // ----------------------

    void SetInsertedString(string insertion) {
        assert(action_ == InsertionSHM);
        insertion_ = insertion;
    }

    string GetInsertedString() const {
        assert(action_ == InsertionSHM);
        return insertion_;
    }

    // ----------------------
    string ToString() const {
        stringstream ss;
        ss << "Mutation at " << position_ << ", type: ";
        if(IsSubstitution())
            ss << "substitution, new nucleotide: " << GetSubstitution();
        else if(IsDeletion())
            ss << "deletion, #deleted: " << GetNumDeleted();
        else if(IsInsertion())
            ss << "insertion, inserted string: " << GetInsertedString();
        return ss.str();
    }
};

ostream& operator<<(ostream& out, const SHM &shm) {
    out << shm.ToString();
    return out;
}

// ----------------------------------------------------------------------------

class SHMSettings {
    vector<SHM> mutations_;
    set<size_t> position_set_;

public:
    void Add(SHM mutation) {
        if(position_set_.find(mutation.Position()) != position_set_.end())
            return;
        mutations_.push_back(mutation);
        position_set_.insert(mutation.Position());
    }

    size_t Size() const { return mutations_.size(); }

    typedef vector<SHM>::iterator shm_iterator;

    shm_iterator begin() { return mutations_.begin(); }

    shm_iterator end() { return mutations_.end(); }

    typedef vector<SHM>::const_iterator shm_citerator;

    shm_citerator cbegin() const { return mutations_.cbegin(); }

    shm_citerator cend() const { return mutations_.cend(); }
};

ostream& operator<<(ostream &out, const SHMSettings &settings) {
    out << "# mutations: " << settings.Size() << endl;
    for(auto it = settings.cbegin(); it != settings.cend(); it++)
        out << *it << endl;
    return out;
}