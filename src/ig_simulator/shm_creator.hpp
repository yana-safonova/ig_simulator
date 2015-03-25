#pragma once

#include "ig_structs/ig_structure_structs.hpp"
#include "repertoire.hpp"

template<class IgVariableRegionPtr, class SHMCreationStrategy, class SHMCreationSettings>
class SHMCreator {
    SHMCreationStrategy strategy_;
public:
    SHMCreator(SHMCreationStrategy strategy) :
        strategy_(strategy) { }

    IgVariableRegionPtr CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMCreationSettings settings =  strategy_.CreateSHM(ig_variable_region_ptr);
        ig_variable_region_ptr->AddSHMSettings(settings);
        return ig_variable_region_ptr;
    }
};

// ----------------------------------------------------------
// Simple SHM strategy

template<class IgVariableRegionPtr, class SHMSettings>
class SHMCreationStrategy {
public:
    SHMSettings CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMSettings settings;
        SHM mutation(0, SubstitutionSHM);
        settings.Add(mutation);
        return settings;
    }
};

typedef SHMCreationStrategy<HC_VariableRegionPtr, SHMSettings> HC_SHMCreationStrategy;
typedef SHMCreationStrategy<LC_VariableRegionPtr, SHMSettings> LC_SHMCreationStrategy;

// ---------------------------------------------------------

typedef SHMCreator<HC_VariableRegionPtr, HC_SHMCreationStrategy, SHMSettings> HC_SHMCreator;
typedef SHMCreator<LC_VariableRegionPtr, LC_SHMCreationStrategy, SHMSettings> LC_SHMCreator;
