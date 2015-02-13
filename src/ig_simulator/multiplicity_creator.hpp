#pragma once

#include "repertoire.hpp"

template<class IgVariableRegionPtr>
class BaseMultiplicityCreator {
public:
    virtual size_t AssignMultiplicity(IgVariableRegionPtr ig_variable_region_ptr) = 0;
};

template<class IgVariableRegionPtr>
class ExponentialMultiplicityCreator : public BaseMultiplicityCreator<IgVariableRegionPtr> {
    double lambda_;

public:
    ExponentialMultiplicityCreator(double lambda) :
            lambda_(lambda) { }

    size_t AssignMultiplicity(IgVariableRegionPtr ig_variable_region_ptr) {
        return size_t((-1) * log(1 - double(rand()) / RAND_MAX) / lambda_) + 1;
    }
};

typedef ExponentialMultiplicityCreator<HC_VariableRegionPtr> HC_ExponentialMultiplicityCreator;
typedef ExponentialMultiplicityCreator<LC_VariableRegionPtr> LC_ExponentialMultiplicityCreator;