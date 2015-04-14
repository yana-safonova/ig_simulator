#pragma once

#include "include_me.hpp"

struct BasicRepertoireParams {
    size_t base_repertoire_size;
    size_t mutated_repertoire_size;
    size_t final_repertoire_size;

    BasicRepertoireParams() :
        base_repertoire_size(0),
        mutated_repertoire_size(0),
        final_repertoire_size(0) { }

    BasicRepertoireParams(size_t base_repertoire_size,
        size_t mutated_repertoire_size,
        size_t final_repertoire_size) :
        base_repertoire_size(base_repertoire_size),
        mutated_repertoire_size(mutated_repertoire_size),
        final_repertoire_size(final_repertoire_size) { }
};

ostream& operator<<(ostream &out, const BasicRepertoireParams &params) {
    out << "Base repertoire size: " << params.base_repertoire_size << endl;
    out << "Expected size of mutated repertoire: " << params.mutated_repertoire_size << endl;
    out << "Expected size of final repertoire: " << params.final_repertoire_size << endl;
    return out;
}

struct PatternSHMParams {
    size_t min_number_pattern_shm;
    size_t max_number_pattern_shm;
    double substitution_propability;

    PatternSHMParams(size_t min_number_pattern_shm,
        size_t max_number_pattern_shm,
        double substitution_propability) :
        min_number_pattern_shm(min_number_pattern_shm),
        max_number_pattern_shm(max_number_pattern_shm),
        substitution_propability(substitution_propability) { }

    PatternSHMParams() :
            min_number_pattern_shm(0),
            max_number_pattern_shm(0),
            substitution_propability(0) { }

    static PatternSHMParams CreateStandardParams() {
        return PatternSHMParams(1, 5, .8);
    }
};

struct CDR_SHMParams {
    size_t min_number_mutations;
    size_t max_number_mutations;
    double mutation_in_fr_prop;

    CDR_SHMParams() :
        min_number_mutations(0),
        max_number_mutations(0),
        mutation_in_fr_prop(0) { }

    CDR_SHMParams(size_t min_number_mutations,
                  size_t max_number_mutations,
                  double mutation_in_fr_prop) :
            min_number_mutations(min_number_mutations),
            max_number_mutations(max_number_mutations),
            mutation_in_fr_prop(mutation_in_fr_prop) { }

    static CDR_SHMParams CreateStandardParams() {
        return CDR_SHMParams(1, 5, .25);
    }
};

struct HC_InputParams {
    string vgenes_fname;
    string dgenes_fname;
    string jgenes_fname;

    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;
    CDR_SHMParams cdr_shm_params;

    HC_InputParams() :
            vgenes_fname(""),
            dgenes_fname(""),
            jgenes_fname(""),
            basic_repertoire_params(),
            pattern_shm_params() { }

    void PrintDatabaseParams() {
        cout << "V gene file: " << vgenes_fname << endl;
        cout << "D gene file: " << dgenes_fname << endl;
        cout << "J gene file: " << jgenes_fname << endl;
    }
};

struct LC_InputParams {
    string vgenes_fname;
    string jgenes_fname;

    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;
    CDR_SHMParams cdr_shm_params;

    LC_InputParams() :
            vgenes_fname(""),
            jgenes_fname(""),
            basic_repertoire_params(),
            pattern_shm_params() { }

    void PrintDatabaseParams() {
        cout << "V gene file: " << vgenes_fname << endl;
        cout << "J gene file: " << jgenes_fname << endl;
    }
};

