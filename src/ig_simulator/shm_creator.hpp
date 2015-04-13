#pragma once

#include "ig_structs/ig_structure_structs.hpp"
#include "repertoire.hpp"
#include "shm_strategies/rgyw_wrcy_strategy.hpp"

template<class IgVariableRegionPtr, class SHMCreationStrategy, class SHMCreationSettings>
class SHMCreator {
    SHMCreationStrategy strategy_;
public:
    SHMCreator(SHMCreationStrategy strategy) :
        strategy_(strategy) { }

    IgVariableRegionPtr CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMCreationSettings settings =  strategy_.CreateSHM(ig_variable_region_ptr);
        cout << "Creation SHM settings" << endl;
        cout << settings;
        ig_variable_region_ptr->SetSHMSettings(settings);
        return ig_variable_region_ptr;
    }
};

// ----------------------------------------------------------
// composite strategy
// todo
// ---------------------------------------------------------


typedef SHMCreator<HC_VariableRegionPtr, HC_RgywWrcySHMStrategy, SHMSettings> HC_SHMCreator;
typedef SHMCreator<LC_VariableRegionPtr, LC_RgywWrcySHMStrategy, SHMSettings> LC_SHMCreator;
