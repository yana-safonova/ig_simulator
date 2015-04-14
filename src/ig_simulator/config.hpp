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
        return PatternSHMParams(2, 5, .8);
    }
};

struct HC_InputParams {
    string vgenes_fname;
    string dgenes_fname;
    string jgenes_fname;

    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;

    HC_InputParams() :
            vgenes_fname(""),
            dgenes_fname(""),
            jgenes_fname(""),
            basic_repertoire_params(),
            pattern_shm_params() { }
};

struct LC_InputParams {
    string vgenes_fname;
    string jgenes_fname;

    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;

    LC_InputParams() :
            vgenes_fname(""),
            jgenes_fname(""),
            basic_repertoire_params(),
            pattern_shm_params() { }
};

