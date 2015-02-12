#pragma once

#include "include_me.hpp"
#include "gene_database.hpp"

class HC_RemovingSettings {
    size_t v_end_len_;
    size_t d_start_len_;
    size_t d_end_len_;
    size_t j_start_len_;

public:
    HC_RemovingSettings():
            v_end_len_(),
            d_start_len_(),
            d_end_len_(),
            j_start_len_() { }

    HC_RemovingSettings(size_t v_end_len, size_t d_start_len, size_t d_end_len, size_t j_start_len) :
        v_end_len_(v_end_len),
        d_start_len_(d_start_len),
        d_end_len_(d_end_len),
        j_start_len_(j_start_len) { }

    size_t VEndLen() const { return v_end_len_; }
    size_t DStartLen() const { return d_start_len_; }
    size_t DEndLen() const { return d_end_len_; }
    size_t JStartLen() const { return j_start_len_; }
};

ostream& operator<<(ostream &out, const HC_RemovingSettings &obj) {
    out << "V end (# removed): " << obj.VEndLen() << endl;
    out << "D start (# removed): " << obj.DStartLen() << endl;
    out << "D end (# removed): " << obj.DEndLen() << endl;
    out << "J start (# removed): " << obj.JStartLen() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class HC_PInsertionSettings {
    size_t v_end_len_;
    size_t d_start_len_;
    size_t d_end_len_;
    size_t j_start_len_;

public:
    HC_PInsertionSettings() :
            v_end_len_(),
            d_start_len_(),
            d_end_len_(),
            j_start_len_() {
    }

    HC_PInsertionSettings(size_t v_end_len, size_t d_start_len, size_t d_end_len, size_t j_start_len) :
        v_end_len_(v_end_len),
        d_start_len_(d_start_len),
        d_end_len_(d_end_len),
        j_start_len_(j_start_len) { }

    size_t VEndLen() const {
        return v_end_len_;
    }

    size_t DStartLen() const {
        return d_start_len_;
    }

    size_t DEndLen() const {
        return d_end_len_;
    }

    size_t JStartLen() const {
        return j_start_len_;
    }
};

ostream& operator<<(ostream &out, const HC_PInsertionSettings &obj) {
    out << "V end (palyndrome length): " << obj.VEndLen() << endl;
    out << "D start (palyndrome length): " << obj.DStartLen() << endl;
    out << "D end (palyndrome length): " << obj.DEndLen() << endl;
    out << "J start (palyndrome length): " << obj.JStartLen() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class HC_NInsertionSettings {
    string vd_insertion_;
    string dj_insertion_;

public:
    HC_NInsertionSettings() :
        vd_insertion_(),
        dj_insertion_() { }

    HC_NInsertionSettings(string vd_insertion, string dj_insertion) :
            vd_insertion_(vd_insertion),
            dj_insertion_(dj_insertion) { }

    string VD_Insertion() const { return vd_insertion_; }

    string DJ_Insertion() const { return dj_insertion_; }
};

ostream& operator<<(ostream &out, const HC_NInsertionSettings &obj) {
    out << "VD insertion: " << obj.VD_Insertion() << ", len: " << obj.VD_Insertion().size() << endl;
    out << "DJ insertion: " << obj.DJ_Insertion() << ", len: " << obj.DJ_Insertion().size() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class HC_VDJ_Recombination {
    const HC_GenesDatabase &db_;
    size_t vgene_index_;
    size_t dgene_index_;
    size_t jgene_index_;

    // endonuclease removing settings
    HC_RemovingSettings removing_settings_;

    // palindromic settings
    HC_PInsertionSettings p_insertion_settings_;

    // n-nucleotides insertion
    HC_NInsertionSettings n_insertion_settings_;

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
        out << "Endonuclease removals:" << endl;
        out << removing_settings_;
        out << "P nucleotides settings:" << endl;
        out << p_insertion_settings_;
        out << "N nucleotides settings:" << endl;
        out << n_insertion_settings_;
    }

    string VgeneName() {
        return db_.GetByIndex(variable_gene, vgene_index_).GeneName();
    }

    string DgeneName() {
        return db_.GetByIndex(diversity_gene, dgene_index_).GeneName();
    }

    string JgeneName() {
        return db_.GetByIndex(join_gene, jgene_index_).GeneName();
    }

    void AddRemovingSettings(HC_RemovingSettings removing_settings) {
        removing_settings_ = removing_settings;
    }

    const HC_RemovingSettings RemovingSettings() const { return removing_settings_; }

    void AddPInsertionSettings(HC_PInsertionSettings p_insertion_settings) {
        p_insertion_settings_ = p_insertion_settings;
    }

    const HC_PInsertionSettings PInsertionSettings() const { return p_insertion_settings_; }

    void AddNInsertionSettings(HC_NInsertionSettings n_insertion_settings) {
       n_insertion_settings_ = n_insertion_settings;
    }

    const HC_NInsertionSettings NInsertionSettings() const { return n_insertion_settings_; }

    const HC_GenesDatabase & GeneDB() const { return db_; }
};

typedef shared_ptr<HC_VDJ_Recombination> HC_VDJ_Recombination_Ptr;

// ----------------------------------------------------------------------------
// todo implement me!

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

typedef shared_ptr<LC_VDJ_Recombination> LC_VDJ_Recombination_Ptr;